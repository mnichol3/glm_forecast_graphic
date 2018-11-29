# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 15:05:28 2018

@author: Salty Pete

Adapted from Alex Chan's code found here
-> https://alexwlchan.net/2017/07/listing-s3-keys/
"""


import boto3
from botocore.handlers import disable_signing
import botocore
import re
import os
import sys

# TODO:
# - Un-hardcode the year & julian day in the prefix var
# - Turn into function & return glm fname list

glm_fnames = []

s3 = boto3.client('s3')
s3.meta.events.register('choose-signer.s3.*', disable_signing) 

keys = []
prefix = 'GLM-L2-LCFA/2018/255/14/OR_GLM-L2-LCFA_G16'
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
    #print(x)
    
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

# return glm_fnames
