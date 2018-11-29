# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 09:20:05 2018

@author: Salty Pete
"""

import pandas as pd
from urllib.request import urlopen
import re
import os

BASE_URL = "https://www.nhc.noaa.gov/archive/recon/"


def hours_mins_2_mins(time):
    if (type(time) != str):
        print("Error: Time must be of type str")
        return -1
    else:
        hrs = int(time[:2])
        mins = int(time[2:])
        
        mins_tot = (hrs * 60) + mins
        
        return mins_tot
        


def dec_2_deg(coord):
    
    splt = coord.split('.')
    
    dec = float(splt[1]) / 60
    deg = float(splt[0]) + dec
    
    if (len(coord) > 5):    # Implies we're dealing with a longitude coordinate
        deg = deg * -1
        
    deg = str(deg)[:6]
    
    return deg



def get_vdm(date, time, storm_name, octant = 'REPNT2'):
    
    vdm_dict = {}
    
    if (type(date) == int):
        print('ERROR: Date argument must be of type str')
        return -1
    
    year = date[:4]
    month = date[4:6]
    day = date[6:]
    
    url = BASE_URL + year + "/" + octant + "/"
    
    try:
        df = pd.read_html(url, skiprows=[0,1,2,3,4,5,6,7,8,9])[0]
    
    except Exception:
        return -1
    
    else:
        df = df.dropna(how="all")
        fnames = df[1].tolist()
        
        obsDateTime = year + month + day
        # Matches all files published during that day for both NOAA & USAF flights
        regex = re.compile(r'(\w{6}-\w{4}\.)' + re.escape(obsDateTime) + r'(\d{4}\.txt)')
        files = list(filter(regex.search, fnames))
        
        matched = False
        
        while (not matched):
            # Get the VDM closest to our given time
            diff = 100000
            min_indx = -1
            int_time = hours_mins_2_mins(time)
            i = 0
            while (i < len(files)):
                curr_diff = abs(int_time - hours_mins_2_mins(files[i][20:24]))
                if (curr_diff < diff):
                    diff = curr_diff
                    min_indx = i
                i += 1
                
            
            # Read the VDM 
            fname = files[min_indx]
            try:
                data = urlopen(url + fname)
                data_bytes = data.read()
                data_str = data_bytes.decode('utf-8')
                data.close()
            
            except Exception as e:
                print("Error fetching file from url in get_vdm_fnames")
                print(e)
                return -1
            
            else:
               # USAF HDBO files have to extra lines at the end containing "$$"
               # and ";". We shell remove these
               data = data_str.splitlines()
               #print(data)
               
               for x in data:
                   if (storm_name.upper() in x):
                       matched = True
               
               # If this is not a VDM for our desired storm, remove the file name
               # from the list and look for a new one
               if (not matched):
                   files.remove(fname)
               else:
                   # Date time format: dd/hh:mm:ss
                   vdm_dict['date_time'] = data[4][3:]
                   coords = data[5].split(' ')
                   vdm_dict['lat'] = coords[1] + coords[3]
                   vdm_dict['lon'] = coords[4] + coords[6]
                   vdm_dict['min SLP'] = coords[7][3:]
                   vdm_dict['inbnd max sfc wind b&r'] = data[12][3:]
                   vdm_dict['inbnd max FL wind b&r'] = data[14][3:]
                   vdm_dict['outbnd max sfc wind b&r'] = data[16][3:]
                   vdm_dict['outbnd max FL wind b&r'] = data[18][3:]
                   vdm_dict['acft storm info'] = data[24][3:]
                   
                   return vdm_dict
                   
                   
    
def vdm_df(date_time_start, date_time_end, storm_name, octant = 'REPNT2'):
    
    if (type(date_time_start) != str):
        date_time_start = str(date_time_start)
        
    if (type(date_time_end) != str):
        date_time_end = str(date_time_end)
    
    date_start = date_time_start[:-4]
    year_start = date_time_start[:4]
    month_start = date_time_start[4:6]
    time_start = date_time_start[-4:]
    date_time_start_int = int(date_time_start)
    
    date_end = date_time_end[:-4]
    month_end = date_time_end[4:6]
    time_end = date_time_end[-4:]
    date_time_end_int = int(date_time_end)
    
    url = BASE_URL + year_start + "/" + octant + "/"
    
    storm_name = storm_name.upper()
    
    try:
        df = pd.read_html(url, skiprows=[0,1,2,3,4,5,6,7,8,9])[0]
    
    except Exception:
        return -1
    
    else:
        df = df.dropna(how="all")
        fnames = df[1].tolist()
        
        obsDateTime = year_start + month_start
        # Matches all files published during that month for both NOAA & USAF flights
        regex1 = re.compile(r'(\w{6}-\w{4}\.)' + re.escape(obsDateTime) + r'(\d{6}\.txt)')
        files1 = list(filter(regex1.search, fnames))
        
        obsDateTime = year_start + month_end
        regex2 = re.compile(r'(\w{6}-\w{4}\.)' + re.escape(obsDateTime) + r'(\d{6}\.txt)')
        files2 = list(filter(regex2.search, fnames))
        
        files = files1 + files2
        #print(files)
        files = [x for x in files if (int(x[12:24]) <= date_time_end_int and 
                                      int(x[12:24]) >= date_time_start_int)]
        #print(files)
        data_list = []
        for x in files:
            url = BASE_URL + year_start + "/" + octant + "/" + x
            curr_list = []
            
            try:
                data = urlopen(url)
                data_bytes = data.read()
                data_str = data_bytes.decode('utf-8')
                data.close()
            
            except Exception as e:
                print("Error fetching file from url in get_vdm_fnames")
                print(e)
                return -1
            
            else:
               print(x)     # Leave this in as a progress indicator
               month = x.split('.')[1][4:6]
               # USAF HDBO files have to extra lines at the end containing "$$"
               # and ";". We shell remove these
               data = data_str.splitlines()
               #print(data)
               
               if (year_start == '2018'):
                   if (storm_name in data[24]):
                       
                       # Date time format: dd/hh:mm:ss
                       date_time = data[4][3:].split('/')       
                       date = month + date_time[0]
                       time = date_time[2].split(':')
                       min_hr = time[0] + time[1]
                       curr_list.append(int(date + min_hr))     # Date time
                       
                       coords = data[5].split(' ')
                       curr_list.append(coords[1] + coords[3])  # lat
                       curr_list.append(coords[4] + coords[6])  # lon
                       mslp = re.search('(\d{3,4})', data[7])
                       if (mslp):
                           curr_list.append(mslp.group(1))      # min SLP
                       else:
                           curr_list.append(0)
                       curr_list.append(data[12][3:])           # inbnd max sfc wind b&r
                       curr_list.append(data[14][3:])           # inbnd max FL wind b&r
                       curr_list.append(data[16][3:])           # outbnd max sfc wind b&r
                       curr_list.append(data[18][3:])           # outbnd max FL wind b&r
                       curr_list.append(data[24][3:])           # acft storm info
                       
                       data_list.append(curr_list)
               else:
                   if (storm_name in data[20]):
                       # Date time format: dd/hh:mm:ss
                       date_time = data[4][3:].split('/')       
                       date = month + date_time[0]
                       time = date_time[1].split(':')
                       min_hr = time[0] + time[1]
                       curr_list.append(int(date + min_hr))     # Date time
                       
                       lats = data[5].split(' ')
                       lat = dec_2_deg(lats[1] + '.' + lats[3])
                       lons = data[6].split(' ')
                       lon = dec_2_deg(lons[2] + '.' + lons[4])
                       curr_list.append(lat)                    # lat
                       curr_list.append(lon)                    # lon
                       mslp = re.search('(\d{3,4})', data[12])
                       if (mslp):
                           curr_list.append(mslp.group(1))      # min SLP
                       else:
                           curr_list.append(0)
                       curr_list.append(data[9][3:])            # inbnd max sfc wind b&r
                       curr_list.append(data[11][3:])           # inbnd max FL wind b&r
                       outbnd = re.search('/\s(\d{1,2})\sNM', data[21])
                       if (outbnd):
                           curr_list.append(outbnd.group(1))    # outbnd max FL wind b&r
                       else:
                           curr_list.append(0)
                       curr_list.append(data[20][3:])           # acft storm info
                       
                       data_list.append(curr_list)
           
        # Dynamically create dataframe column names
        alpha = 'ABCDEFGHIJ'
        col_names = []
        i = 0
        while (i < len(data_list[0])):
            col_names.append(alpha[i])
            i += 1
        #col_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        vdm_df = pd.DataFrame(data_list, columns=col_names)
        vdm_df = vdm_df.sort_values(by=['A'])
        
        new_fname = "VDM-" + storm_name + "-" + date_time_start + "-" + date_time_end + ".txt"

        if not os.path.exists(os.path.join(octant, year_start + '-' + storm_name)):
            os.makedirs(os.path.join(octant, year_start + '-' + storm_name))
                   
        path = os.path.join(octant, year_start + '-' + storm_name, new_fname)
        with open(path, 'w') as file:
            vdm_df.to_csv(file, sep = ',', header=False, index=False)

        return vdm_df



def read_vdm_csv(fname):
    
    if ('2018' in fname):
        col_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    else:
        col_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        
    vdm_df = pd.read_csv(fname, sep=",", header=None)
    vdm_df.columns = col_names
    
    return vdm_df



def main():
    date = '20180913'
    time = '1600'
    storm_name = 'florence'
    #vdm_dict = get_vdm(date, time, storm_name)
    #print(vdm_dict)
    vdm_df('201708300000', '201709130600', 'irma')
    
    
    
    
if __name__ == "__main__":
    main() 