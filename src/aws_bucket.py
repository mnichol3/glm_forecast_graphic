# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 15:05:28 2018

@author: Matt Nicholson

Adapted from Alex Chan's code found here
-> https://alexwlchan.net/2017/07/listing-s3-keys/
"""


import boto3
from botocore.handlers import disable_signing
import botocore
import re
import os
import sys



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



def glm_dl(date):
    """
    Download GOES-16 GLM data files from NOAA's Amazon AWS server

    Parameters
    ------------
    date : str
        Date & time of the desired files, in a 1-hr block. Format: YYYYMMDDHH

    Returns
    ------------
    glm_fnames : list of str
        List of downloaded GLM filenames as they are stored locally
        Format: YYYMMDDHHMM.nc
    """
    glm_fnames = []

    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    year = date[:4]
    hour = date[-2:]
    julian_day = calc_julian_day(date)

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

    if (get_os() == 'linux'):
        path = '/home/mnichol3/Documents/senior-rsch/data/glm'
    else:
        path = 'D:\Documents\senior-research-data\glm'
    dl_count = 0

    for x in keys:

        fname_match = re.search('s(\d{14})', x)

        if (fname_match):
            fname = fname_match[1] + '.nc'
            try:
                sys.stdout.write("\rDownloading GLM datafile: {}".format(fname))
                s3.download_file('noaa-goes16', x, os.path.join(path, fname))
                glm_fnames.append(fname)
                dl_count += 1
                sys.stdout.flush()
            except botocore.exceptions.ClientError as e:
                if e.response['Error']['Code'] == "404":
                    print("The file does not exist in this AWS bucket.")
                else:
                    print("Error downloading file from AWS")

    sys.stdout.write("\rFinished! Files downloaded: {}              ".format(dl_count))
    sys.stdout.flush()

    return glm_fnames



def abi_dl(date, sector):
    """
    Downlad GOES-16 ABI data files from NOAA's Amazon AWS server

    Parameters
    ------------
    date : str
        Date & time of the desired files, in a 1-hr block. Format: YYYYMMDDHH
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

    if (get_os() == 'linux'):
        path = '/home/mnichol3/Documents/senior-rsch/data/abi'
    else:
        path = 'D:\Documents\senior-research-data\abi'
    dl_count = 0

    for x in keys:

        fname_match = re.search('s(\d{11})', x)

        if (fname_match):
            fname = fname_match[1] + '.nc'
            try:
                sys.stdout.write("\rDownloading GLM datafile: {}".format(fname))
                s3.download_file('noaa-goes16', x, os.path.join(path, fname))
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

    return abi_fnames


if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <aws_bucket> as main...')
