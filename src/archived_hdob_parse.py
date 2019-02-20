# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 18:52:14 2018

@author: Matt Nicholson

These functions identify, download, & parse High Density Observation (HDOB) files
on the National Hurricane Center's data server

Example usage:

date = '20181009'
time = '1458'
storm_name = 'michael'

file_dict = get_archived_fnames(date)
#write_fname_dict(file_dict, date, storm_name, 1)
#print(file_dict['usaf'])
nfile_dict = create_archived_HDOB_file(file_dict['usaf'], storm_name, date+time)
print(nfile_dict)
#write_fname_dict(nfile_dict, date, storm_name, 2)
plot_flight(nfile_dict)
#print(create_archived_HDOB_file(['TEST'], 'GORDON'))
"""

import pandas as pd
from urllib.request import urlopen
from pandas.compat import StringIO
import os
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
from . import manager as mg
#import manager as mg

BASE_URL = 'https://www.nhc.noaa.gov/archive/recon/'
PATH_LINUX = '/home/mnichol3/Documents/senior-rsch/data/hdob'
PATH_WIN = r'D:\Documents\senior-research-data\hdob'

def minutes_degrees(coord, kywrd):
    """
    Converts coordinates from decimal minutes to decimal degrees

    Parameters
    ------------
    coord : str
        A single coordinate, either longitude or latitude, with a single char
        at the end indicating the hemisphere
    kywrd : str
        Indicates whether the coordinate is a latitude or longitude coordinate,
        as they are formatted differently

    Returns
    ------------
    deg : str
        Decimal degree coordinate
    """
    if (kywrd == 'lat'):
        deg = coord[:2]
        mins = coord[2:4]
    elif (kywrd == 'lon'):
        deg = coord[:3]
        mins = coord[3:5]
    else:
        print('Coordinate Transform error: Invalid keyword')
        return -1
    dec = float(mins) / 60
    deg = float(deg) + dec

    if (coord[-1] == "S" or coord[-1] == "W"):
        deg = deg * -1

    deg = str(deg)[:6]

    return deg



def create_archived_HDOB_file(fnames, storm_name_param, date_time):
    """
    This function takes a list of HDOB filenames and downloads them from the
    NHC aircraft recon data server. It also creates a text file holding all the
    processed data in a subdirectory named after the octant in which the obs
    were taken. Directory structure will be similar to that of the NHC server

    Ex:
    dir_name = octant
    new_fname =  octant + "-" + obs_date + "-" + acft_callsign + "-" + storm_name + ".txt"

    Parameters
    ------------
    fnames : list of strs
        List of strings of files to download
    storm_name_param : str
        Name of the storm the aircraft obs were taken in, used only for file
        naming purposes
    date_time : str
        Date & time of desired files. Format: YYYYMMDDHH ?

    Returns
    ------------
    new_fnames : list of str
        List of filenames of processed HDOB files, stored locally
    """
    global BASE_URL
    new_fnames = {}

    if (mg.get_os() == 'linux'):
        path = PATH_LINUX
    else:
        path = PATH_WIN

    # NHC data server url format:
    # https://www.nhc.noaa.gov/archive/recon/2018/AHONT1/AHONT1-KNHC.201801032351.txt

    file_date_time = [int(x.split('.')[1]) for x in fnames]
    #print(file_date_time)

    nrst_time = min(file_date_time, key=lambda x:abs(x-int(date_time)))

    nrst_time_indx = file_date_time.index(nrst_time)

    # First we'll move backwards
    i = nrst_time_indx
    skipped = 0
    while (i >= 0 and skipped < 5):
        fname = fnames[i]


    #for fname in fnames:

        date = fname[12:20]
        year = date[:4]

        octant = fname[:6]

        # TODO: UNCOMMENT URL
        url = BASE_URL + year + "/" + octant + "/" + fname
        print(url)
        #url = "https://www.nhc.noaa.gov/archive/recon/2018/AHONT1/AHONT1-KWBC.201809091621.txt"
        #url = "https://www.nhc.noaa.gov/archive/recon/2018/AHONT1/AHONT1-KNHC.201809041753.txt"
        #print('\n!!! WARNING: updateHDOB url is hardcoded !!!\n')

        try:
            data = urlopen(url)
            data_bytes = data.read()
            data_str = data_bytes.decode('utf-8')
            data.close()

        except Exception as e:
            print("Error fetching file from url in create_archived_HDOB_file")
            print(e)
            return -1

        else:
           # USAF HDBO files have to extra lines at the end containing "$$"
           # and ";". We shell remove these
           if (len(data_str.splitlines()) > 23):
                # Remove last three lines
                data = data_str.splitlines()[:-3]
           else:
                data = data_str.splitlines()

           file_date = data[2].split(' ')[2][:2]
           obs_time = data[2].split(' ')[2][2:]

           line_strings = [x for x in data[3].split(' ') if x != '']

           acft_callsign = line_strings[0]
           storm_name = line_strings[2]
           obs_date = line_strings[-1]
           obs_num = line_strings[-2]

           obs_date = obs_date[:6] + file_date

           # Set flags for first iteration
           if (i == nrst_time_indx):
               prev_calls = acft_callsign
               prev_obs_num = 45

           if (storm_name.lower() == storm_name_param.lower() and (prev_calls == acft_callsign)
               and (int(obs_num) < int(prev_obs_num))):

               names = ["ObsTime", "Lat", "Lon", "StatAirPress", "GeoPotHgt",
                        "SfcPressDVal", "T_air", "T_dew", "WndDirSpd", "WndPk",
                        "SfcWndPk", "RainRate", "Qflags"]

               df = pd.read_csv(StringIO('\n'.join(data)), delimiter=r"\s+",
                                dtype = str, names = names, skiprows=[0, 1, 2, 3])

               df['Lat'] = df.apply(lambda row: minutes_degrees(row['Lat'], 'lat'), axis=1)
               df['Lon'] = df.apply(lambda row: minutes_degrees(row['Lon'], 'lon'), axis=1)

               new_fname = storm_name + "-" + obs_date + "-" + obs_time + "-" + acft_callsign + "-" + obs_num + ".txt"

               prev_calls = acft_callsign
               prev_obs_num = obs_num

               if not os.path.exists(os.path.join(path, octant, year + '-' + storm_name)):
                   os.makedirs(os.path.join(path, octant, year + '-' + storm_name))

               abs_path = os.path.join(path, octant, year + '-' + storm_name, new_fname)
               with open(abs_path, 'w') as file:
                   df.to_csv(file, sep = ',', header=False, index=False)

               new_key = obs_date + '-' + storm_name + '-' + acft_callsign
               if (new_key in new_fnames.keys()):
                   new_fnames[new_key].append(new_fname)
               else:
                   new_fnames[new_key] = [new_fname]
           else:
                skipped += 1

        i -= 1

    # Now iterate forwards in time
    i = nrst_time_indx
    skipped = 0
    while (i < len(fnames) and skipped < 5):
        fname = fnames[i]


    #for fname in fnames:

        date = fname[12:20]
        year = date[:4]

        octant = fname[:6]

        # TODO: UNCOMMENT URL
        url = BASE_URL + year + "/" + octant + "/" + fname
        print(url)
        #url = "https://www.nhc.noaa.gov/archive/recon/2018/AHONT1/AHONT1-KWBC.201809091621.txt"
        #url = "https://www.nhc.noaa.gov/archive/recon/2018/AHONT1/AHONT1-KNHC.201809041753.txt"
        #print('\n!!! WARNING: updateHDOB url is hardcoded !!!\n')

        try:
            data = urlopen(url)
            data_bytes = data.read()
            data_str = data_bytes.decode('utf-8')
            data.close()

        except Exception as e:
            print("Error fetching file from url in create_archived_HDOB_file")
            print(e)
            return -1

        else:
           # USAF HDBO files have to extra lines at the end containing "$$"
           # and ";". We shell remove these
           if (len(data_str.splitlines()) > 23):
                # Remove last three lines
                data = data_str.splitlines()[:-3]
           else:
                data = data_str.splitlines()

           file_date = data[2].split(' ')[2][:2]
           obs_time = data[2].split(' ')[2][2:]

           line_strings = [x for x in data[3].split(' ') if x != '']

           acft_callsign = line_strings[0]
           storm_name = line_strings[2]
           obs_date = line_strings[-1]
           obs_num = line_strings[-2]

           obs_date = obs_date[:6] + file_date

           # Set flags for first iteration
           if (i == nrst_time_indx):
               prev_calls = acft_callsign
               prev_obs_num = 0

           if (storm_name.lower() == storm_name_param.lower() and (prev_calls == acft_callsign)
               and (int(obs_num) > int(prev_obs_num))):

               names = ["ObsTime", "Lat", "Lon", "StatAirPress", "GeoPotHgt",
                        "SfcPressDVal", "T_air", "T_dew", "WndDirSpd", "WndPk",
                        "SfcWndPk", "RainRate", "Qflags"]

               df = pd.read_csv(StringIO('\n'.join(data)), delimiter=r"\s+",
                                dtype = str, names = names, skiprows=[0, 1, 2, 3])

               df['Lat'] = df.apply(lambda row: minutes_degrees(row['Lat'], 'lat'), axis=1)
               df['Lon'] = df.apply(lambda row: minutes_degrees(row['Lon'], 'lon'), axis=1)

               new_fname = storm_name + "-" + obs_date + "-" + obs_time + "-" + acft_callsign + "-" + obs_num + ".txt"

               prev_calls = acft_callsign
               prev_obs_num = obs_num

               if not os.path.exists(os.path.join(path, octant, year + '-' + storm_name)):
                   os.makedirs(os.path.join(path, octant, year + '-' + storm_name))

               abs_path = os.path.join(path, octant, year + '-' + storm_name, new_fname)
               with open(abs_path, 'w') as file:
                   df.to_csv(file, sep = ',', header=False, index=False)

               new_key = obs_date + '-' + storm_name + '-' + acft_callsign
               if (new_key in new_fnames.keys()):
                   new_fnames[new_key].append(new_fname)
               else:
                   new_fnames[new_key] = [new_fname]
           else:
                skipped += 1

        i += 1

    return new_fnames



def get_archived_fnames(date, octant = 'AHONT1'):
    """
    Fetches NOAA & USAF aircraft observation files published on the given date,
    as well as the day before & day after in order to capture the full length
    of flights the occured in whole or in part on the specified date

    Parameters
    ------------
    date : str
        Desired date of the aircraft obs. Format: YYYYMMDD
    octant : str
        Octant code indicating where the observation was taken. Default is
        Atlantic Basin

    Returns
    ------------
    fname_dict : dictionary; keys : str, values : list of str
        A dictionary containing a list of NOAA aircraft observations and a list
        of USAF aircraft observtions. EX: {'usaf': usaf_list, 'noaa': noaa_list}
    """
    global BASE_URL

    # filename format: AHONT1-KWBC.201809091451.txt
    #                              YYYYMMDDzzzz

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

        obsDateTime = year + month
        # Matches all files published during that month for both NOAA & USAF flights
        regex = re.compile(r'(\w{6}-\w{4}\.)' + re.escape(obsDateTime) + r'(\d{6}\.txt)')

        # Adjust the date to the day before
        # Doesnt handle leap year case
        # or new year case
        # But if theres recon flights on new years we have bigger problems...
        prev_day = int(day) - 1
        prev_mon = month

        if (prev_day < 1):
            prev_mon = int(prev_mon)
            prev_mon -= 1

            if (prev_mon in [1, 3, 5, 7, 8, 10, 12]):
                prev_day = 31
            elif (prev_mon in [4, 6, 9, 11]):
                prev_day = 30
            else:
                prev_day = 28

            if (prev_mon < 10):
                prev_mon = '0' + str(prev_mon)
            else:
                prev_mon = str(prev_mon)
            prev_day = str(prev_day)
        else:
            if (prev_day < 10):
               prev_day = '0' + str(prev_day)
            else:
                prev_day = str(prev_day)


        prev_date = year + prev_mon + prev_day

        # Adjust date of following day
        next_day = int(day) + 1
        next_month = month

        if (next_day >= 30):
            next_month = int(next_month)

            if (next_month in [4, 6, 9, 11] and next_day > 30):
                next_month += 1
                next_day = 1
            elif (next_month in [1, 3, 5, 7, 8, 10, 12] and next_day > 31):
                next_month += 1
                next_day = 1

            # New years case. Not sure why a flight would occur on NYE
            if (next_month > 12):
                year = int(year) + 1
                year = str(year)
                next_month = 1

            if (next_month < 10):
                next_month = '0' + str(next_month)
            else:
                next_month = str(next_month)

        if (next_day < 10):
            next_day = '0' + str(next_day)
        else:
            next_day = str(next_day)

        next_date = year + next_month + next_day


        # What this ends up doing is presenting the last 2 days of files
        files = list(filter(regex.search, fnames))

        # Get all files for the current date, plus files published after (before)
        # 12z on the previous (next) date to ensure the entire duration of a flight
        # that occurs on the given date is oncluded
        files = [x for x in files if (date in x) or ((prev_date + '1') in x) or
                 ((next_date + '0') in x)]

        if (len(files) > 0):

            usaf = [x for x in files if 'KNHC' in x]

            noaa = [x for x in files if 'KWBC' in x]


#            print("USAF:")
#            print(usaf)
#            print("\n")
#            print("NOAA:")
#            print(noaa)
        fname_dict = {'usaf': usaf, 'noaa': noaa}

        return fname_dict



def plot_flight(fname_dict):
    """
    Plots HDOB wind speeds & directions as color-coded wind barbs for a given
    flight

    Parameters
    ------------
    fname_dict : dictionary; keys : str, values : list of str
        Dictionary of something

    Returns
    ------------
    Plot of flight's HDOBs, with wind speeds & directions as color-coded wind barbs
    """

    df_list = []

    if (mg.get_os() == 'linux'):
        path = PATH_LINUX
    else:
        path = PATH_WIN

    for key, val in fname_dict.items():

        #print("Key: " + key)
        #print(val)
        info = val[0].split('-')
        storm_name = info[0]
        obs_date = info[1]
        year = obs_date[:4]
        callsign = info[3]

        base_path = os.path.join(path, "AHONT1", year + "-" + storm_name)

        for x in val:
            abs_path = os.path.join(base_path, x)
            curr = pd.read_csv(abs_path, sep=",", dtype = str, header=None)
            curr.columns = ["ObsTime", "Lat", "Lon", "StatAirPress", "GeoPotHgt",
                            "SfcPressDVal", "T_air", "T_dew", "WndDirSpd", "WndPk",
                            "SfcWndPk", "RainRate", "Qflags"]

            df_list.append(curr)

        data = pd.concat(df_list, ignore_index=True)
        #print(data)

        lats = data['Lat'].astype(np.float32).tolist()
        lons = data['Lon'].astype(np.float32).tolist()

        winds = data['WndDirSpd'].tolist()

        #        wind dir       wind spd
        winds = [(int(x[:3]), int(x[3:])) for x in winds if '///' not in x]

        u = [x[1] * math.sin(math.radians(x[0])) for x in winds]
        v = [x[1] * math.cos(math.radians(x[0])) for x in winds]

        wind_spd = [x[1] for x in winds]

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor((0, 0, 0))
        ax.coastlines(resolution='10m', color='gray')

        cmap = cm.tab20c
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
        bounds = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        m = cm.ScalarMappable(cmap=cmap, norm=norm)
        m.set_array(wind_spd)

        ax.barbs(lons, lats, u, v, wind_spd, length=6, pivot='middle',
                 sizes=dict(emptybarb=0.25, spacing=0.2, height=0.5),
                 linewidth=0.6, cmap=cmap, norm=norm)

        x_min, x_max, y_min, y_max = plt.axis('square')

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray',
                          alpha=0.5, linestyle='--', draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right=False

        xticks = np.linspace(x_min, x_max, 6)
        xticks = np.concatenate([[x_min], xticks, [x_max]])
        yticks = np.linspace(y_min, y_max, 6)
        yticks = np.concatenate([[y_min], yticks, [y_max]])
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)

        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
        gl.ylabel_style = {'color': 'red', 'weight': 'bold'}

        cbar = plt.colorbar(m, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds, pad=0.02)
        cbar.ax.set_ylabel('Flight-level Wind Speed (kts)')

        initial_date_time = val[0].split('-')
        idt = initial_date_time[1] + ' ' + initial_date_time[2] + 'z'

        final_date_time = val[-1].split('-')
        fdt = final_date_time[1] + ' ' + final_date_time[2] + 'z'

        plt.title('Aircraft Obs - ' + callsign, loc='left', fontweight='bold')
        plt.title(storm_name.upper() + ' - ' + idt + ' to ' + fdt, loc='right')

        plt.tight_layout()

        plt.show()

if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <archived_hdob_parse> as main...')
