# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 12:26:21 2019

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



def padding_zero(int_to_pad, lim):
    if (int_to_pad < lim):
        result = '0' + str(int_to_pad)
    else:
        result = str(int_to_pad)
    return result



###############################################################################
# Splits a start date_time & end date_time into 1-hr chunks
#
# @param    start_date_time (str)       Beginning of the date_time chunk.
#                                       Format: YYYYMMDDHH or YYYYMMDDHHMM
#
# @param    end_date_time (str)         End of the date_time chunk.
#                                       Format: YYYYMMDDHH or YYYYMMDDHHMM
#
# @return   chunks (list of str)        List of date_times between start_date_time
#                                       & end_date_time (inclusive), split into
#                                       1-hr chunks 
# Ex -> print(date_time_chunk('2018091218', '2018091305'))
###############################################################################
def date_time_chunk(start_date_time, end_date_time):
    
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



###############################################################################
# Calculates the julian day (days since 01 Jan) from the date param. Needed to 
# download from AWS server
#
# @param    date (str)          Date to convert to julian day. Format: YYYYMMDD
#                               HH and/or MM (minute) may be included at the end but not
#                               required (e.g., YYYYMMDDHHMM)
#
# @return   julian_day (str)    Converted julian day date. 
###############################################################################
def calc_julian_day(date):
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



###############################################################################
# Downloads GOES-16 GLM data files from NOAA's AWS server
#
# @param    date (str)                  Date & time of the desired files. Time 
#                                       is 1-hr block. Format: YYYYMMDDHH
#
# @return   glm_fnames (list of str)    List of filenames downloaded. 
#                                       Format: YYYYMMDDHHMM.nc
###############################################################################
def glm_dl(date):
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



###############################################################################
# Downloads GOES-16 ABI data files from NOAA's AWS server
#
# @param    date (str)                  Date & time of the desired files. Time 
#                                       is 1-hr block. Format: YYYYMMDDHH
#
# @param    sector (str)                M1 -> meso 1, M2 -> meso 2, C -> conus
#
# @return   abi_fnames (list of str)    List of filenames downloaded. 
#                                       Format: YYYYMMDDHHMM.nc
###############################################################################
def abi_dl(date, sector):
    
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

