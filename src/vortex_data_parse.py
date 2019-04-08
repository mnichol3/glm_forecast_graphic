# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 09:20:05 2018

@author: Matt Nicholson

These functions deal with downloading, parsing, & processing NOAA & USAF
Vortex Data Message (VDM) aircraft observations

"""

import pandas as pd
from urllib.request import urlopen
import re
import os
import sys
import numpy as np
#from . common import get_os, padding_zero, date_time_chunk
from common import get_os, padding_zero, date_time_chunk
import config


BASE_URL = 'https://www.nhc.noaa.gov/archive/recon/'
PATH_LINUX_VDM = '/media/mnichol3/easystore/data/vdm'
PATH_WIN = r'D:\Documents\senior-research-data\vdm'



def calc_min_list(date_time_start, date_time_end, mode='time_list'):
    """
    Returns a list of datetime strings for every minute between date_time_start
    & date_time_end (inclusive)
    Date time format: MoMoDDHHMM (Mo = month, M = minute)

    Parameters
    ----------
    date_time_start : str
        Beginning date & time. Must be before date_time_end
        Format : MoMoDDHHMM (Mo = month, M = minute)

    date_time_end : str
        Ending date & time. Must be after date_time_start
        Format : MoMoDDHHMM (Mo = month, M = minute)

    mode : str
        String indicating what to return.

        mode = 'time_list' --> A list of datetime strings occuring between
                               date_time_start & date_time_end (inclusive) will
                               be returned, in 1-minute increments

        mode = 'tot_mins' --> An int representing the total minutes from
                              date_time_start to date_time_end (inclusive)
                              will be returned


    Returns
    -------
    times : list of str
        List containing the datetime strings occuring between date_time_start
        & date_time_end (inclusive)

    OR

    tot_mins : int
        An int representing the total minutes from date_time_start to date_time_end
        (inclusive). Counting starts at 1 (date_time_start = minute 1)

    Notes
    -----
    Could just return len(times) instead of tot_mins, but that would take up
    a lot of memory for no reason

    """
                    # J   F   M   A   M   J   J   A   S   O   N   D
                    # 1   2   3   4   5   6   7   8   9   10  11  12
                    # 0   1   2   3   4   5   6   7   8   9   10  11
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    times = []
    tot_mins = 0

    curr = date_time_start

    while (curr != date_time_end):

        # Not very efficient but saves memory
        if (mode == 'time_list'):
            times.append(curr)
        elif (mode == 'tot_mins'):
            tot_mins += 1
        else:
            print('ERROR: Invalid mode argument (vortex_data_parse.calc_min_list)')
            sys.exit(0)

        month = curr[:2]
        day = curr[2:4]
        hr = curr[4:6]
        mins = int(curr[6:8])

        mins += 1

        if (mins > 59):
            mins = 0
            hr = int(hr)
            hr += 1

            if (hr > 23):
                hr = 0
                day = int(day)
                day += 1

                if (day > days_per_month[int(month) - 1]):
                    day = 1
                    month = int(month)
                    month += 1
                    month = str(month)
                    month = padding_zero(month, 10)
                day = padding_zero(day, 10)
            hr = padding_zero(hr, 10)
        mins = padding_zero(mins, 10)

        curr = month + day + hr + mins

    if (mode == 'time_list'):
        times.append(date_time_end)
        return times
    else:
        tot_mins += 1
        return tot_mins



def hours_mins_2_mins(time):
    """
    Converts a time consisting of hours & minutes to minutes

    Parameters
    ------------
    time : str
        Time, in hours & minutes, to be converted to minutes


    Returns
    ------------
    mins_tot : int
        Time converted from hours:minutes to minutes
    """
    if (type(time) != str):
        print("Error: Time must be of type str")
        return -1
    else:
        hrs = int(time[:2])
        mins = int(time[2:])

        mins_tot = (hrs * 60) + mins

        return mins_tot



def dec_2_deg(coord):
    """
    Converts a coordinate from decimal minutes to decimal degrees

    Parameters
    ------------
    coord : str
        Coordinate to be converted from decimal minutes to decimal degrees


    Returns
    ------------
    deg : str
        Coordinate transformed into decimal degrees
    """

    splt = coord.split('.')

    dec = float(splt[1]) / 60
    deg = float(splt[0]) + dec

    if (len(coord) > 5):    # Implies we're dealing with a longitude coordinate
        deg_str = '-' + str(deg)
        return deg_str[:7]
    else:
        return str(deg)[:6]

    #return deg



def get_single_vdm(date, time, storm_name, octant = 'REPNT2'):
    """
    Downloads a VDM from within a specific storm for a given day & time.

    Parameters
    ------------
    date : str
        Date the observation was taken
        FORMAT: YYYYMMDD

    time : str
        Time the observation was taken
        FORMAT: YYYYMMDD

    storm_name : str
        Name of the storm the observation was taken in

    octant : str
        Octant that the obs was taken within. Default is REPNT2


    Returns
    ------------
    vdm_dict : dictionary of VDM data

    Ex:
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

    Notes
    -----
    * ONLY WORKS FOR 2018 VDM COORDINATE FORMAT
    """

    vdm_dict = {}
    year = date[:4]
    month = date[4:6]
    day = date[6:]
    url = BASE_URL + year + "/" + octant + "/"

    if (type(date) == int):
        print('ERROR: Date argument must be of type str')
        return -1

    if (int(year) < 2018):
        print('ERROR: get_vdm function only works for years after 2017')
        print('get_vdm may fail on 2018 VDMs published before June also')
        sys.exit(0)

    try:
        df = pd.read_html(url, skiprows=[0,1,2,3,4,5,6,7,8,9], index_col=0)[0]

    except Exception as e:
        print(e)
        sys.exit(0)

    else:
        fnames = [f for f in df.iloc[:,0].tolist() if str(f) != 'nan']
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
                sys.exit(0)

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

                   curr_data = data[6].split(' ')
                   vdm_dict['min SLP'] = curr_data[1] + ' mb'
                   vdm_dict['inbnd max sfc wind b&r'] = data[12][3:]
                   vdm_dict['inbnd max FL wind b&r'] = data[14][3:]
                   vdm_dict['outbnd max sfc wind b&r'] = data[16][3:]
                   vdm_dict['outbnd max FL wind b&r'] = data[18][3:]
                   vdm_dict['acft storm info'] = data[24][3:]

                   return vdm_dict



def vdm_df(date_time_start, date_time_end, storm_name, octant = 'REPNT2'):
    """
    Creates a Pandas dataframe with accumulated VDM data, writes the dataframe
    to a csv, and returns the dataframe

    Parameters
    ------------
    date_time_start : str
        Starting date & time of VDMs
        FORMAT: YYYYMMDDHHMM

    date_time_end : str
        Ending date & time of VDMs
        FORMAT: YYYYMMDDHHMM

    storm_name : str
        Name of the storm the VDMs were taken in

    octant : str
        Octant the VDMs were taken in


    Returns
    ------------
    vdm_df : Pandas dataframe
        A Pandas dataframe containing an accumulation of data from multiple VDMs
        The dataframe is then written to a csv.

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

    Notes
    -----
    * USAF RPENT2 VDM's published after May 2018 have coordinates in decimal degrees
    * USAF RPENT2 VDM's published in 2017 have coordinates in decimal minutes
    * NOAA RPENT2 VDM's published in 2017 have coordinates in decimal minutes

    """

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
        df = pd.read_html(url, skiprows=[0,1,2,3,4,5,6,7,8,9], index_col=0)[0]

    except Exception:
        print("Error reading url in vdm_df")
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
            print('Oops! No vortex data messages published between ' + date_time_start
                    + ' and ' + date_time_start)
            sys.exit(0)

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
                sys.exit(0)

            else:
               print("Downloading: " + x)     # Leave this in as a progress indicator
               month = x.split('.')[1][4:6]
               # USAF HDBO files have to extra lines at the end containing "$$"
               # and ";". We shall remove these
               data = data_str.splitlines()
               #print(data)

               if (year_start == '2018'):
                   if (storm_name in data[24]):

                       # Date time format: dd/hh:mm:ss
                       date_time = data[4][3:].split('/')
                       date = month + date_time[0]
                       time = date_time[1].split(':')
                       min_hr = time[0] + time[1]
                       curr_list.append(date + min_hr)     # Date time

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
                       curr_list.append(date + min_hr)     # Date time

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
        vdm_df = pd.DataFrame(data_list, columns=col_names, dtype=str)
        vdm_df = vdm_df.sort_values(by=['A'])
        vdm_df = vdm_df.drop_duplicates(subset=['A', 'B', 'C'], keep='last')

        new_fname = "VDM-" + storm_name + "-" + date_time_start + "-" + date_time_end + ".txt"

        if (get_os() == 'linux'):
            path = PATH_LINUX_VDM
        else:
            path = PATH_WIN

        if not os.path.exists(os.path.join(path, octant, year_start + '-' + storm_name)):
            os.makedirs(os.path.join(path, octant, year_start + '-' + storm_name))

        abs_path = os.path.join(path, octant, year_start + '-' + storm_name, new_fname)
        with open(abs_path, 'w') as file:
            vdm_df.to_csv(file, sep = ',', header=False, index=False)

        return vdm_df



def read_vdm_csv(fname):
    """
    Reads a csv file that a Pandas dataframe of accumulated VDM data was written
    to & returns the csv data as a Pandas dataframe

    Parameters
    ------------
    fname : str
        Name of the csv file the Pandas dataframe was written to


    Returns
    ------------
    vdm_df : Pandas dataframe
        Dataframe containing the data in the fname csv.

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
    """

    if ('2018' in fname):
        col_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    else:
        col_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    vdm_df = pd.read_csv(fname, sep=",", header=None, dtype=str)
    vdm_df.columns = col_names
    vdm_df = vdm_df.drop_duplicates(subset=['A'], keep='last')

    return vdm_df



def track_interp(vdm_df, year, interval):
    """
    This function takes a dataframe of accumulated VDM data and linearly
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
    delta_ts = []
    rmw_vals = []

    obs_times = vdm_df['A'].tolist()

    lats = [float(i[:-1]) for i in vdm_df['B'].tolist()]
    lons = [float(j[:-1]) for j in vdm_df['C'].tolist()]

    min_slp = [float(k) for k in vdm_df['D'].tolist()]

    rmw_in = vdm_df['F'].tolist()
    rmw_in = [x.split(' ')[2] for x in rmw_in]

    if (year == '2018'):
        rmw_out = vdm_df['H'].tolist()
    else:
        rmw_out = vdm_df['G'].tolist()

    rmw_out = [x.split(' ')[2] for x in rmw_out]

    if (len(rmw_in) != len(rmw_out)):
        print("Error parsing inbound & outbound RMW in vdp.track_interp")
        sys.exit(0)
    else:
        lim = len(rmw_in)
        idx = 0
        while (idx < lim):
            #avg_rmw = (int(rmw_in[idx]) + int(rmw_out[idx])) / 2
            #rmw_vals.append(avg_rmw)
            rmw_vals.append(max(int(rmw_in[idx]), int(rmw_out[idx])))
            idx += 1

    if (interval == "hour"):

        t_0 = obs_times[1]
        t_f = obs_times[-1]

        for t in obs_times[1:]:
            delta_ts.append(calc_min_list(t_0, t, 'tot_mins'))

        lats = lats[1:]
        lons = lons[1:]

        date_times = date_time_chunk(year + t_0[:-2], year + t_f[:-2])
        date_times = date_times[1:]
        time_span = np.zeros(len(date_times))

        # Create the times (in minutes from t_0) that we want interpolated
        # coordinates for
        for i, val in enumerate(date_times):
            time_span[i] = calc_min_list(t_0, val[4:], 'tot_mins')

        lat_interp = np.interp(time_span, delta_ts, lats)
        lon_interp = np.interp(time_span, delta_ts, lons)
        rmw_interp = np.interp(time_span, delta_ts, rmw_vals[1:])

        #min_slp_interp = np.interp(time_span, delta_ts, min_slp[1:])

        """
        ######### For debugging purposes #########
        print("Inc xp: " + np.all(np.diff(time_span) > 0))
        print('time_span len: ' + str(len(time_span)))
        print('date_times len: ' + str(len(date_times)))
        print('lons len: ' + str(len(lons)))
        print('lats len: ' + str(len(lats)))
        for x in zip(lon_interp, lat_interp):
            print(x)
        ##########################################
        """


    elif (interval == "min"):

        t_0 = obs_times[0]
        t_f = obs_times[-1]

        # Determine the time between t_0 and each of the other obs, in minutes
        for t in obs_times:
            delta_ts.append(calc_min_list(t_0, t, 'tot_mins'))

        time_span = np.arange(0, delta_ts[-1])

        lat_interp = np.interp(time_span, delta_ts, lats)
        lon_interp = np.interp(time_span, delta_ts, lons)
        rmw_interp = np.interp(time_span, delta_ts, rmw_vals)

        #min_slp_interp = np.interp(time_span, delta_ts, min_slp)

        date_times = calc_min_list(t_0, t_f, 'time_list')

        """
        ######### For debugging purposes #########
        print("Inc xp: " + np.all(np.diff(time_span) > 0))
        print('time_span len: ' + str(len(time_span)))
        print('date_times len: ' + str(len(date_times)))
        print('lons len: ' + str(len(lons)))
        print('lats len: ' + str(len(lats)))
        ##########################################
        """



    else:
        print('ERROR: Invalid interval parameter (vortex_data_parse.track_interp)')
        sys.exit(0)


    processed_vdm = pd.DataFrame({'date_time': date_times,
                                  'time_span': time_span,
                                  'lons': lon_interp,
                                  'lats': lat_interp,
                                  'rmw' : rmw_interp
    })

    return processed_vdm



if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <vortex_data_parse> as main...')
