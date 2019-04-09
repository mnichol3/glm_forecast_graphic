"""
Created on Mon Apr 8 09:20:05 2018

@author: Matt Nicholson

These functions deal with downloading, parsing, & processing NOAA & USAF
TEMP DROP (dropsonde) aircraft observations

"""

import pandas as pd
from urllib.request import urlopen
import re
import os.makedirs
from os.path import join, isfile, exists
import sys
import numpy as np
#from . common import get_os, padding_zero, date_time_chunk
from common import get_os, padding_zero, date_time_chunk
import config



BASE_URL = 'https://www.nhc.noaa.gov/archive/recon/'
PATH_LINUX_TMPDRP = '/media/mnichol3/easystore/data/tmpdrp'


def tmpdrp_df(date_time_start, date_time_end, storm_name, octant = 'REPNT3'):
    """
    Creates a Pandas dataframe with accumulated dropsonde data, writes the dataframe
    to a csv, and returns the dataframe

    Parameters
    ------------
    date_time_start : str
        Starting date & time of dropsondes
        FORMAT: YYYYMMDDHHMM

    date_time_end : str
        Ending date & time of dropsondes
        FORMAT: YYYYMMDDHHMM

    storm_name : str
        Name of the storm the dropsondes were taken in

    octant : str
        Octant the dropsondes were taken in


    Returns
    ------------
    tmpdrp_df : Pandas dataframe
        A Pandas dataframe containing an accumulation of data from multiple
        dropsondes. The dataframe is then written to a csv.

        Columns:

            datetime    : Date & time of the observation. Format: MoMoDDHHMM
            lat         : Latitude at which the dropsonde was released. Format: XX.X
            lon         : Longitude at which the dropsonde was released. Format: XXX.X
                          *Note: Will include leading '-' in W Hemisphere (default octant)
            dlm_wind    : Deep Layer Mean Wind, in knots. Format: ddfff
                            dd:  True direction from which the wind is blowing,
                                 rounded to the nearest 5 degrees. Reported in
                                 hundreds and tens digits. The unit digit (0 or 5)
                                 is added to the hundreds digit of wind speed
                            fff: Wind speed in knots. Hundreds digit is sum of
                                 speed and unit digit of direction, i.e. 295 deg
                                 at 125 knots is encoded as 29625
    """
    data_list = []

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

    if (get_os() != 'linux'):
        print("Dude, Windows sucks")
        sys.exit(0)

    new_fname = "tmpdrp-" + storm_name + "-" + date_time_start + "-" + date_time_end + ".txt"

    abs_path = join(PATH_LINUX_TMPDRP, octant, year_start + '-' + storm_name, new_fname)

    # See if accumulated dropsonde data file already exists. If it does, just read
    # it and return it rather than re-download everything
    if isfile(abs_path):
        tmpdrp_df = read_tmpdrp_csv(abs_path)

    else:

        try:
            df = pd.read_html(url, skiprows=[0,1,2,3,4,5,6,7,8,9], index_col=0)[0]

        except Exception:
            print("Error reading url (temp_drop_parse.tmpdrp_df)")
            print(e)
            sys.exit(0)

        else:
            fnames = [f for f in df.iloc[:,0].tolist() if str(f) != 'nan']

            obsDateTime = year_start + month_start
            # Matches all files published during that month for both NOAA & USAF flights
            regex1 = re.compile(r'(\w{6}-\w{4}\.)' + re.escape(obsDateTime) + r'(\d{6}\.txt)')
            files1 = list(filter(regex1.search, fnames))

            obsDateTime = year_start + month_end
            regex2 = re.compile(r'(\w{6}-\w{4}\.)' + re.escape(obsDateTime) + r'(\d{6}\.txt)')
            files2 = list(filter(regex2.search, fnames))

            files = files1 + files2

            files = [x for x in files if (int(x[12:24]) <= date_time_end_int and
                                          int(x[12:24]) >= date_time_start_int)]

            if (not files):
                print('Oops! No dropsonde data published between ' + date_time_start
                        + ' and ' + date_time_start)
                sys.exit(0)

            for x in files:
                url = BASE_URL + year_start + "/" + octant + "/" + x
                curr_list = []

                try:
                    data = urlopen(url)
                    data_bytes = data.read()
                    data_str = data_bytes.decode('utf-8')
                    data.close()

                except Exception as e:
                    print("Error fetching file from url (temp_drop_parse.tmpdrp_df)")
                    print(e)
                    sys.exit(0)

                else:
                   print("Downloading: " + x)    # Leave this in as a progress indicator
                   month = x.split('.')[1][4:6]
                   # USAF HDBO files have to extra lines at the end containing "$$"
                   # and ";". We shall remove these

                   data = data_str.splitlines()

                   if (storm_name in data_str):

                       print("Processing: " + x + '\n')

                       # Idx = 2
                       date_time = data[2].split(' ')[2] # DDHHMM
                       date = month + date_time[:2]
                       time = date_time[2:]
                       curr_list.append(date + time)     # MoMoDDHHMM

                       # Idx = 3
                       curr_line = data[3].split(' ')
                       curr_line = [x for x in curr_line if x != '']

                       # Latitude
                       if (curr_line[2][:2] == '99'):
                           curr_list.append(curr_line[2][2:4] + '.' + curr_line[2][4:])
                       else:
                           print('Error parsing dropsonde latitude in ' + x + '
                                 (temp_drop_parse.tmpdrp_df)')
                           sys.exit(0)

                       # Longitude
                       lon = float(curr_line[3][1:4] + '.' + curr_line[3][4:]) * -1
                       curr_list.append(str(lon))


                       dlm_wnd = re.search('DLM WND (\d{5})', data_str).group(1)

                       if (dlm_wnd):
                           curr_list.append(dlm_wnd.group(1))
                       else:
                           print("Error parsing DLM WND (temp_drop_parse.tmpdrp_df)")
                           sys.exit(0)

                       data_list.append(curr_list)

            col_names = ['datetime', 'lat', 'lon', 'dlm_wind']
            tmpdrp_df = pd.DataFrame(data_list, columns=col_names, dtype=str)
            tmpdrp_df = tmpdrp_df.sort_values(by=['A'])
            tmpdrp_df = tmpdrp_df.drop_duplicates(subset=col_names, keep='last')

            #new_fname = "tmpdrp-" + storm_name + "-" + date_time_start + "-" + date_time_end + ".txt"

            file_path = join(PATH_LINUX_TMPDRP, octant, year_start + '-' + storm_name)
            if not exists(file_path):
                makedirs(file_path)

            #abs_path = join(file_path, new_fname) # Defined at top of function
            with open(abs_path, 'w') as file:
                tmpdrp_df.to_csv(file, sep = ',', header=False, index=False)

    return tmpdrp_df



def read_tmpdrp_csv(fname):
    """
    Reads a csv file that a Pandas dataframe of accumulated VDM data was written
    to & returns the csv data as a Pandas dataframe

    Parameters
    ------------
    fname : str
        Name of the csv file the Pandas dataframe was written to

    Returns
    ------------
    tmpdrp_df : Pandas dataframe
        A Pandas dataframe containing an accumulation of data from multiple
        dropsondes.

        Columns:

            datetime    : Date & time of the observation. Format: MoMoDDHHMM
            lat         : Latitude at which the dropsonde was released. Format: XX.X
            lon         : Longitude at which the dropsonde was released. Format: XXX.X
                          *Note: Will include leading '-' in W Hemisphere (default octant)
            dlm_wind    : Deep Layer Mean Wind, in knots. Format: ddfff
                            dd:  True direction from which the wind is blowing,
                                 rounded to the nearest 5 degrees. Reported in
                                 hundreds and tens digits. The unit digit (0 or 5)
                                 is added to the hundreds digit of wind speed
                            fff: Wind speed in knots. Hundreds digit is sum of
                                 speed and unit digit of direction, i.e. 295 deg
                                 at 125 knots is encoded as 29625
    """
    col_names = ['datetime', 'lat', 'lon', 'dlm_wind']

    tmpdrp_df = pd.read_csv(fname, sep=",", header=None, dtype=str)
    tmpdrp_df.columns = col_names
    tmpdrp_df = tmpdrp_df.drop_duplicates(subset=col_names, keep='last')

    return tmpdrp_df



def tmpdrp_interp(tmpdrp_df, datetimes=None):
    """
    This function takes a dataframe of accumulated dropsonde data and linearly
    interpolates the time and low pressure center coordinates on 1-minute
    intervals. The times used are measured as minutes since the first published
    observation for that storm

    Parameters
    ----------
    vdm_df : Pandas DataFrame
        DataFrame containing accumulated VDM data

        Columns:

            2018+ files:
                A : Date time
                B : Latitude
                C : Longitude
                D : Minimum SLP
                E : Inbound maximum surface wind, incl. bearing & range
                F : Inbound maximum flight-level wind, incl. bearing & range
                G : Outbound maximum surface wind, incl. bearing & range
                H : Outbound maximum flight-level wind, incl. bearing & range
                I : Aircfraft info

            Files previous to 2018:
                A : Date time
                B : Latitude
                C : Longitude
                D : Minimum SLP
                E : Inbound maximum surface wind, incl. bearing & range
                F : Inbound maximum flight-level wind, incl. bearing & range
                G : Outbound maximum flight-level wind, incl. bearing & range
                H : Aircfraft info

    year : str
        Year the observations were taken

    interval : str
        String indicating whether to interpolate to 1-minute or 1-hour


    Returns
    -------
    interp_vdf : Pandas Dataframe
        Pandas Dataframe holding the interpolated data


    Notes
    -----


    """




tmpdrp_df('201809010900', '201809140300', 'FLORENCE')
