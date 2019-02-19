# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 12:26:21 2019

@author: Matt Nicholson

This file contains scripts used to download GOES-16 Geostationary Lightning
Mapper (GLM) & Advanced Baseline Imager (ABI) data files from NOAA's Amazon
Web Service (AWS) bucket.

Adapted from Alex Chan's code found here
-> https://alexwlchan.net/2017/07/listing-s3-keys/
"""


import boto3
from botocore.handlers import disable_signing
import botocore
import re
import os
import sys
import csv



def get_os():
    """
    Determines the type of operating system being used. Needed for when we are
    loading & saving local files later

    Parameters
    ------------
    none

    Returns
    ------------
    os_type : str
        Type of OS the script is running on
        Ex: 'linux' is the script is running on a Linux OS, such as Ubuntu
    """
    os_type = sys.platform
    return os_type



def padding_zero(int_to_pad, lim):
    """
    Adds a leading, or padding, zero to an integer if needed. Returns the int
    as a str

    Parameters
    ------------
    int_to_pad : int
        The integer to add padding/leading zero to, if needed

    lim : int
        Used to check to see if the param int needs a padding zero added to it.
        Ex:
        int_to_pad = 11, lim = 10 -> no padding needed, '11' is returned
        int_to_pad = 9, lim = 10 -> padding needed, '09' is returned

    Returns
    ------------
    result : str
        int_to_pad cast as a string, regardless if a padding 0 was added to it
    """
    if (int_to_pad < lim):
        result = '0' + str(int_to_pad)
    else:
        result = str(int_to_pad)
    return result



def date_time_chunk(start_date_time, end_date_time):
    """
    Splits a start date_time & end date_time into 1-hr chunks

    Ex:
    date_time_chunk('2018091218', '2018091305')

    yields:
    ['2018091218', '2018091219', '2018091220', ... , '2018091304', '2018091305']

    Parameters
    ------------
    start_date_time : str
        Beginning of the date_time chunk. Format: YYYYMMDDHH or YYYYMMDDHHMM

    end_date_time : str
        End of the date_time chunk. Format: YYYYMMDDHH or YYYYMMDDHHMM

    Returns
    ------------
    chunks : list of str
        List of date_times between start_date_time & end_date_time (inclusive),
        split into 1-hr chunks
    """

                    # J   F   M   A   M   J   J   A   S   O   N   D
                    # 1   2   3   4   5   6   7   8   9   10  11  12
                    # 0   1   2   3   4   5   6   7   8   9   10  11
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    leap_years = [2016, 2020, 2024, 2028, 2032]

    chunks = [start_date_time]

    curr_date_time = start_date_time

    while (curr_date_time != end_date_time):
        curr_year = int(curr_date_time[:4])
        curr_month = int(curr_date_time[4:6])
        curr_day = int(curr_date_time[6:8])
        curr_hour = int(curr_date_time[8:10])

        curr_hour += 1
        if (curr_hour > 23):
            curr_hour = 0
            curr_day += 1
            if (curr_day > 30):
                if (days_per_month[curr_month - 1] == 30):
                    curr_day = 0
                    curr_month += 1
            elif (curr_day > 28 and curr_month == 2):
                if (curr_year not in leap_years):
                    curr_day = 0
                    curr_month += 1
                elif (curr_day > 29):
                    curr_day = 0
                    curr_month += 1

        if (curr_month > 12):
            curr_month = 0
            curr_year += 1

        curr_year = str(curr_year)

        curr_month = padding_zero(curr_month, 10)
        curr_day = padding_zero(curr_day, 10)
        curr_hour = padding_zero(curr_hour, 10)

        curr_date_time = curr_year + curr_month + curr_day + curr_hour

        chunks.append(curr_date_time)

    return chunks



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
    given date, time, & hour.

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
        path = '/home/mnichol3/Documents/senior-rsch/data/glm'
    else:
        path = 'D:\Documents\senior-research-data\glm'

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
                break

        for x in keys:

            fname_match = re.search('s(\d{14})', x)

            if (fname_match):
                fname = fname_match[1] + '.nc'
                local_fname = year + month + day + hour + fname[-8:-3] + '.nc'
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



def abi_dl(date_lst, sector):
    """
    Downlad GOES-16 ABI data files from NOAA's Amazon AWS server

    Parameters
    ------------
    date_lst : list of str
        List of the dates & times of the desired files, in a 1-hr block.
        Format: YYYYMMDDHH

    sector : str
        Sector of the GOES-16 imagery to download. "M1" -> mesoscale 1,
        "M2" -> mesoscale 2, "C" -> CONUS

    Returns
    ------------
    abi_fnames : list of str
        List of filenames downloaded as they are stored locally.
        Format: YYYYMMDDHHMM.nc
    """
    abi_fnames = []
    dl_count = 0

    if (get_os() == 'linux'):
        path = '/home/mnichol3/Documents/senior-rsch/data/abi'
    else:
        path = 'D:\Documents\senior-research-data\glm'

    # Make sure we have a list
    if (type(date_lst) == str):
        date_lst = [date_lst]

    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    keys = []
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
        return
    for date in date_lst:
        year = date[:4]
        hour = date[-2:]
        julian_day = calc_julian_day(date)

        prefix = 'ABI-L1b-RadM/' + year + '/' + julian_day + '/' + hour + '/OR_ABI-L1b-Rad' + sector_prefix + '-M3C01_G16'
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
                break

        for x in keys:

            fname_match = re.search('s(\d{11})', x)

            if (fname_match):
                fname = fname_match[1] + '.nc'
                local_fname = year + month + day + hour + fname[-8:-3] + '.nc'
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
