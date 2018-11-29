# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 13:10:14 2018

@author: Salty Pete
"""

import re

'''
a = 'MAX FL WIND 102 KT 305 / 25 NM 16:44:30Z'

outbnd = re.search('/\s(\d{1,2})\sNM', 'MAX FL WIND 102 KT 305 / 25 NM 16:44:30Z')

if (outbnd):
    print('yes')
'''
    
'''
test_str = '_s20182551459400_'

fname = re.match('_s(\d{14})_', test_str)

if (fname):
    print(fname[1])
'''

import sys
import time

count = 0

while(count < 10):
    sys.stdout.write("\rDownloading GLM datafile: {}".format(count))
    sys.stdout.flush()
    time.sleep(0.5)
    count += 1
sys.stdout.write("\rFinished!                             ")
sys.stdout.flush()

