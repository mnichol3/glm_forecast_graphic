# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 12:26:21 2019

@author: Matt Nicholson

This file contains scripts used to download GOES-16 Geostationary Lightning
Mapper (GLM) & Advanced Baseline Imager (ABI) data files from NOAA's Amazon
Web Service (AWS) bucket.

Code that deals with AWS buckt adapted from Alex Chan's code found here
-> https://alexwlchan.net/2017/07/listing-s3-keys/

--------------------------------------------------------------------------------
NOAA AWS ABI file format:

OR_ABI-L2-CMIPC-M3_G16_sYYYYJJJHHMMSST_eYYYYJJJHHMMSST_cYYYYJJJHHMMSST.nc

Where:
    OR      : Indicates the system is operational
    ABI     : Instrument type (Advanced Baseline Imager, in this case)
    L2      : Indicated Level 2 data
    CMIP    : Cloud and Moisture Imagery Product
    C       : Indicates CONUS file
    M3      : Scan mode
    G16     : Indicates satellite (GOES-16 in this case)
    s#####  : Scan start date & time
    e#####  : Scane end date & time
    c#####  : File creation date & time

    ** Filename date & time formats:
    YYYY : Year
    JJJ  : Julian day
    HH   : Hour, in 23-hr format with leading zero if < 10
    MM   : Minutes, with leading zero if < 10
    SS   : Seconds, with leading zero if < 10
    T    : Tenth second, no leading zero or decimal

Ex:

OR_ABI-L2-MCMIPC-M3_G16_s20181781922189_e20181781924562_c20181781925075.nc

Scan start      : 2018 178 19:22:18:9
Scan end        : 2018 178 19:24:56:2
File created    : 2018 178 19:25:07:5

--------------------------------------------------------------------------------

NOAA AWS GLM Filename format:

<sensor>-<level>-<product short name>/<year>/<julian day>/<hour>/OR_<sensor>
-<level>-<product short name>-M<scanning mode>-C<channel>-G<GOES Satellite>
-s<start time>_e<end time>_c<central time>.nc

"""


import boto3
from botocore.handlers import disable_signing
import botocore
import re
import os
import sys
import csv
from common import get_os, padding_zero, date_time_chunk



PATH_LINUX_ABI = '/home/mnichol3/Documents/senior-rsch/data/abi'
PATH_LINUX_GLM = '/home/mnichol3/Documents/senior-rsch/data/glm'
PATH_WIN_ABI = r'D:\Documents\senior-research-data\abi'
PATH_WIN_GLM = r'D:\Documents\senior-research-data\glm'



def dec_min(date_time):
    """
    Decreases the minute of date_time and adjusts the hour and julian day
    if needed

    Parameters
    ----------
    date_time : str
        Date_time to be adjusted.
        Format: YYYYJdJdJdHHMM, where:
            YYYY   : Year
            JdJdJd : 3-digit julian day
            HH     : Hours
            MM     : Minutes

    Returns
    -------
    adj_dt : str
        Adjusted date_time string

    """
    month = date_time[:2]
    day = date_time[4:7]
    hr = date_time[7:9]
    mins = int(date_time[-2:])

    mins -= 1

    if (mins < 0):
        mins = 0
        hr = int(hr)
        hr -= 1

        if (hr < 0):
            hr = 0
            day = int(day)
            day -= 1

            day = padding_zero(day, 10)
        hr = padding_zero(hr, 10)
    mins = padding_zero(mins, 10)

    adj_dt = month + day + hr + mins

    return adj_dt



def calc_julian_day(date):
    """
    Calculates the Julian Day (days since 01 Jan) for a given date. Needed to
    download files from AWS server

    Parameters
    ------------
    date : str
        Date to be converted to julian day. Format: YYYYMMDD, YYYYMMDDHH, or
        YYYYMMDDHHMM

    Returns
    ------------
    julian_day : str
        Date converted to Julian Day
    """
    # Format: YYYYMMDDHHMM

                    # J   F   M   A   M   J   J   A   S   O   N   D
                    # 1   2   3   4   5   6   7   8   9   10  11  12
                    # 0   1   2   3   4   5   6   7   8   9   10  11
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    leap_years = [2016, 2020, 2024, 2028, 2032]

    year = date[:4]
    month = date[4:6]
    day = date[6:8]

    if (int(year) in leap_years):
        days_per_month[1] += 1

    curr_month = int(month) - 2
    julian_day = int(day)

    while (curr_month >= 0):
        julian_day += days_per_month[curr_month]
        curr_month -= 1

    julian_day = str(julian_day)

    if (int(julian_day) < 100):
        julian_day = '0' + julian_day

        if (int(julian_day) < 10):
            julian_day = '0' + julian_day

    return julian_day



def glm_dl(date_lst):
    """
    Download ALL GOES-16 GLM data files from NOAA's Amazon AWS server for the
    given date  & hour.

    NOAA AWS GLM Filename format:

    <sensor>-<level>-<product short name>/<year>/<julian day>/<hour>/OR_<sensor>
    -<level>-<product short name>-M<scanning mode>-C<channel>-G<GOES Satellite>
    -s<start time>_e<end time>_c<central time>.nc

    13 Feb 19 : Now writes fnames to glm_fnames.csv

    Ex:
    date = '2018092812' (12z Sept 9, 2018)
    All GLM files from 201809281200 (1200z) to 201809281259 (1259z) will be
    downloaded, totaling 180 files as GLM files are published every 20 seconds

    Parameters
    ------------
    date_lst : list of str
        List of date & time string of the desired files, in a 1-hr block.
        Format: YYYYMMDDHH

    Returns
    ------------
    glm_fnames : list of str
        List of downloaded GLM filenames as they are stored locally
        Format: YYYMMDDHHMM.nc
    """
    glm_fnames = []
    dl_count = 0

    if (get_os() == 'linux'):
        path = PATH_LINUX_GLM
    else:
        path = PATH_WIN_GLM

    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    # If a string is passed instead of a list of strings
    if (type(date_lst) == str):
        date_lst = [date_lst]

    for date in date_lst:

        year = date[:4]
        hour = date[-2:]
        julian_day = calc_julian_day(date)
        month = date[4:6]
        day = date[6:8]

        keys = []
        prefix = 'GLM-L2-LCFA/' + year + '/' + julian_day + '/' + hour + '/OR_GLM-L2-LCFA_G16'
        suffix = ''
        kwargs = {'Bucket': 'noaa-goes16', 'Prefix': prefix}

        while True:
            resp = s3.list_objects_v2(**kwargs)
            for obj in resp['Contents']:
                key = obj['Key']
                if key.endswith(suffix):
                    keys.append(key)

            try:
                kwargs['ContinuationToken'] = resp['NextContinuationToken']
            except KeyError:
                break # LEAVE AS A BREAK !!!!!!!

        for x in keys:

            fname_match = re.search('e(\d{14})', x) # Matches scan end date time

            if (fname_match):
                # Filename with month & day as julian day
                fname = fname_match[1] + '.nc'
                # Filename with month & day as month & day
                local_fname = year + month + day + fname[-8:-3] + '.nc'
                try:
                    sys.stdout.write("\rDownloading GLM datafile: {}".format(local_fname))
                    s3.download_file('noaa-goes16', x, os.path.join(path, local_fname))
                    glm_fnames.append(fname)
                    dl_count += 1
                    sys.stdout.flush()
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        print("The file does not exist in this AWS bucket.")
                    else:
                        print("Error downloading file from AWS\n")

    sys.stdout.write("\rFinished! Files downloaded: {}              ".format(dl_count))
    sys.stdout.flush()
    print("\n")
    '''
    fname = 'glm_fnames-' + date_lst[0] + '-' + date_lst[-1] + '.csv'
    with open(os.path.join(path, fname), 'w') as f:
        writer = csv.writer(f, lineterminator='\n')
        for x in glm_fnames:
            writer.writerow([x])
    '''
    return glm_fnames



def abi_dl_multi(date_lst, sector, band=1):
    """
    Download multiple GOES-16 ABI L2 data files from NOAA's Amazon AWS server for
    a given date, hour, and sector

    Parameters
    ------------
    date_lst : list of str
        List of the dates & times of the desired files, in a 1-hr block.
        Format: YYYYMMDDHH

    sector : str
        Sector of the GOES-16 imagery to download. "M1" -> mesoscale 1,
        "M2" -> mesoscale 2, "C" -> CONUS

    band : int
        Integer representing the band to download. Default is 01 (Red visible)
        Range: 1-16


    Returns
    ------------
    abi_fnames : list of str
        List of filenames downloaded as they are stored locally.
        Format: YYYYMMDDHHMM.nc
    """
    abi_fnames = []
    dl_count = 0
    band = padding_zero(band, 10)

    if (get_os() == 'linux'):
        path = PATH_LINUX_ABI
    else:
        path = PATH_WIN_ABI

    # Make sure we have a list
    if (type(date_lst) == str):
        date_lst = [date_lst]

    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    # RadM1 for Meso sector 1
    # RadM2 for Meso sector 2
    if (sector == 'meso1'):
        sector_prefix = 'M1'
    elif (sector == 'meso2'):
        sector_prefix = 'M2'
    elif (sector == 'conus'):
        sector_prefix = 'C'
    else:
        print('Error: Invalid sector parameter!')
        sys.exit(0)

    for date in date_lst:
        year = date[:4]
        hour = date[-2:]
        julian_day = calc_julian_day(date)

        keys = []
        prefix = ('ABI-L2-CMIP' + sector_prefix[0] + '/' + year + '/' + julian_day + '/' + hour +
                  '/OR_ABI-L2-CMIP' + sector_prefix + '-M3C' + band + '_G16')
        suffix = ''
        kwargs = {'Bucket': 'noaa-goes16', 'Prefix': prefix}

        while True:
            resp = s3.list_objects_v2(**kwargs)
            for obj in resp['Contents']:
                key = obj['Key']
                if key.endswith(suffix):
                    keys.append(key)

            try:
                kwargs['ContinuationToken'] = resp['NextContinuationToken']
            except KeyError:
                break # LEAVE AS A BREAK !!!!!!!

        for x in keys:

            fname_match = re.search('e(\d{11})', x) # Matches scan end date time

            if (fname_match):
                fname = fname_match[1] + '.nc'
                local_fname = year + month + day + hour + fname[-7:-3] + '.nc'
                try:
                    sys.stdout.write("\rDownloading ABI datafile: {}".format(fname))
                    s3.download_file('noaa-goes16', x, os.path.join(path, local_fname))
                    abi_fnames.append(fname)
                    dl_count += 1
                    sys.stdout.flush()
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        print("The file does not exist in this AWS bucket.")
                    else:
                        print("Error downloading file from AWS")
                    sys.exit(0)

    sys.stdout.write("\rFinished! Files downloaded: {}              ".format(dl_count))
    sys.stdout.flush()
    print("\n")
    '''
    fname = 'abi_fnames-' + date_lst[0] + '-' + date_lst[-1] + '.csv'
    with open(os.path.join(path, fname), 'w') as f:
        writer = csv.writer(f, lineterminator='\n')
        for x in abi_fnames:
            writer.writerow([x])
    '''
    return abi_fnames



def abi_dl(date_time, sector, band=1):
    """
    Downloads a single GOES-16 ABI L2 data file from NOAA's Amazon AWS server for
    a given datetime & sector

    Parameters
    ------------
    date_time : str
        Date time of the desired file
        Format: YYYYMMDDHHMM

    sector : str
        Sector of the GOES-16 imagery to download. "M1" -> mesoscale 1,
        "M2" -> mesoscale 2, "C" -> CONUS

    band : int
        Integer representing the band to download. Default is 01 (Red visible)
        Range: 1-16


    Returns
    ------------
    abi_fname : str
        Filename downloaded as it is stored locally.
        Format: YYYYMMDDHHMM.nc
    """
    abi_fnames = []
    keys = []
    dl_count = 0
    tries = 0
    success = False
    band = padding_zero(band, 10)

    if (get_os() == 'linux'):
        path = PATH_LINUX_ABI
    else:
        path = PATH_WIN_ABI

    # Make sure we have a list
    if (len(date_time) != 12 ):
        print('ERROR: Invalid date-time string (aws_dl.abi_dl)')
        sys.exit(0)

    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    # RadM1 for Meso sector 1
    # RadM2 for Meso sector 2
    if (sector == 'meso1'):
        sector_prefix = 'M1'
    elif (sector == 'meso2'):
        sector_prefix = 'M2'
    elif (sector == 'conus'):
        sector_prefix = 'C'
    else:
        print('ERROR: Invalid sector parameter (aws_dl.abi_dl)')
        sys.exit(0)


    year = date_time[:4]
    month = date_time[4:6]
    day = date_time[6:8]
    hour = date_time[8:10]
    mins = date_time[-2:]
    julian_day = calc_julian_day(date_time)

    prefix = ('ABI-L2-CMIP' + sector_prefix[0] + '/' + year + '/' + julian_day + '/' + hour +
              '/OR_ABI-L2-CMIP' + sector_prefix + '-M3C' + band + '_G16')
    suffix = ''
    kwargs = {'Bucket': 'noaa-goes16', 'Prefix': prefix}

    while True:
        resp = s3.list_objects_v2(**kwargs)
        for obj in resp['Contents']:
            key = obj['Key']
            if key.endswith(suffix):
                keys.append(key)

        try:
            kwargs['ContinuationToken'] = resp['NextContinuationToken']
        except KeyError:
            break # LEAVE AS A BREAK !!!!!!!

    while (not success and tries < 120):
        for x in keys:

            # Matches scan end date-time
            fname_match = re.search('e(' + year + julian_day + hour + mins + ')', x)

            if (fname_match):
                # eYYYYJJJHHMMSST
                fname = fname_match[1] + '.nc'
                local_fname = year + month + day + hour + fname[-7:-3] + '.nc'
                try:
                    print("Downloading ABI datafile: " + fname)
                    s3.download_file('noaa-goes16', x, os.path.join(path, local_fname))
                    success = True
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        print("The file does not exist in this AWS bucket.")
                    else:
                        print("Error downloading file from AWS")
                    sys.exit(0)

        if (not success):
            print("Oops! No ABI file for that date and time found")
            print("Adjusting and trying again")

            tries += 1

            adj_dt = dec_min(year + julian_day + hour + mins)

            year = adj_dt[:4]
            julian_day = adj_dt[4:7]
            hour = adj_dt[7:9]
            mins = adj_dt[-2:]


    if (success):
        print("Finished!")
        return local_fname
    else:
        print("No ABI files near that date and time were found!")
        return -1



if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <aws_dl> as main...')
