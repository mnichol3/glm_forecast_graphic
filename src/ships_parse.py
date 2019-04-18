# -*- coding: utf-8 -*-
"""
Created on Wed 10 Apr 2019 16:16

@author: Matt Nicholson

These functions obtain and parse Statistical Hurricane Intensity Predictions
Scheme (SHIPS) files


SHIPS filename format:
    YYMMDDHHBBSSYY_ships.txt
    Where:
        YY : Year
        MM : Month
        DD : Day
        HH : Hour
        BB : Basin (AL = Atlantic, EP = Eastern Pacific)
        SS : Storm number
        YY : Year

    ex:
        18102712AL1618_ships.txt

        27 Oct 2018 12z AL16

"""
import urllib.request
import sys
import pandas as pd
from os.path import join, exists, isfile
from os import makedirs, listdir
import csv



BASE_URL = 'ftp://ftp.nhc.noaa.gov/atcf/stext'
PATH_LINUX_SHIPS = '/media/mnichol3/easystore/data/SHIPS'



def fetch_file_local(datetime, storm_name, basin='AL', write=False):
    """
    Locates and parses a SHIPS file for the given datetime & storm

    Parameters
    ----------
    datetime : str
        Date & time of the desired SHIPS data. Format: YYYYMoMoDDHHMM
    storm_name : str
        Name of the storm the desired SHIPS file is related to
    basin : str
        Basin in which the storm is located. AL = Atlantic, EP = Eastern Pacific


    Returns
    -------
    data_dict : dictionary of type key : str, value : str
        Dictionary containing relevant SHIPS data. Keys: datetime, shear_spd,
        shear_dir, sst

    """
    html = None
    data = None
    hour_diff = 0
    data_dict = {}
    synop_times = ['00', '06', '12', '18']
    storm_nums = {'FLORENCE': '06', 'MICHAEL': '14'}
    base_path = PATH_LINUX_SHIPS + "/ftp"

    storm_name = storm_name.upper()
    if (storm_name not in list(storm_nums.keys())):
        print('Error: Storm name not in storm number dictionary (ships_parse.fetch_files)')
        sys.exit(0)

    fnames = [f for f in listdir(base_path) if isfile(join(base_path, f))]

    if (fnames):

        storm_name = storm_name.upper()

        year = datetime[2:4]
        month = datetime[4:6]
        day = datetime[6:8]
        hour_true = datetime[8:10]

        # Determine the hour for the filename
        if (hour_true not in synop_times):
            if (hour_true < '18'):
                prev_idx = next(idx for idx, value in enumerate(synop_times) if value > hour_true)
                prev_idx -= 1
                hour_adj = synop_times[prev_idx]
            else:
                hour_adj = '18'

            hour_diff = int(hour_true) - int(hour_adj)
            hour = hour_adj
        else:
            hour = hour_true

        target_fname = year + month + day + hour + basin + storm_nums[storm_name] + year + '_ships.txt'

        ftp_fname = join(base_path, target_fname)

        print("Processing SHIPS file: " + ftp_fname + "\n")

        with open(ftp_fname, 'r') as response:
            data = response.read()

        if (data):
            data = data.splitlines()

            idx = 2
            while (idx < 14):
                line = data[idx].split(' ')
                line = [val for val in line if val != '']

                if (idx == 2):
                    if (storm_name in line):
                        data_dict['datetime'] = datetime
                    else:
                        print("Error: Storm names do not match (ships_parse.fetch_files)")
                        sys.exit(0)
                elif (idx == 10):
                    if (hour_diff == 0):
                        data_dict['shear_spd'] = line[2]
                    else:
                        data_dict['shear_spd'] = line[3]
                elif (idx == 12):
                    if (hour_diff == 0):
                        data_dict['shear_dir'] = line[2]
                    else:
                        data_dict['shear_dir'] = line[3]
                elif (idx == 13):
                    if (hour_diff == 0):
                        data_dict['sst'] = line[2]
                    else:
                        data_dict['sst'] = line[3]
                idx += 1
        else:
            print("Error retrieving SHIPS file data (ships_parse.fetch_files)")
            sys.exit(0)


    else:
        print("Error retrieving SHIPS filenames (ships_parse.fetch_files)")
        sys.exit(0)

    return data_dict



def fetch_file(datetime, storm_name, basin='AL', write=False):
    """
    Locates and parses a SHIPS file for the given datetime & storm

    Parameters
    ----------
    datetime : str
        Date & time of the desired SHIPS data. Format: YYYYMoMoDDHHMM
    storm_name : str
        Name of the storm the desired SHIPS file is related to
    basin : str
        Basin in which the storm is located. AL = Atlantic, EP = Eastern Pacific
    write : bool
        If true, writes the dataframe to a csv. Default: False


    Returns
    -------
    data_dict : dictionary of type key : str, value : str
        Dictionary containing relevant SHIPS data. Keys: datetime, shear_spd,
        shear_dir, sst

    """
    html = None
    data = None
    hour_diff = 0
    data_dict = {}
    synop_times = ['00', '06', '12', '18']
    storm_nums = {'FLORENCE': '06', 'MICHAEL': '14'}

    storm_name = storm_name.upper()
    if (storm_name not in list(storm_nums.keys())):
        print('Error: Storm name not in storm number dictionary (ships_parse.fetch_files)')
        sys.exit(0)


    with urllib.request.urlopen(BASE_URL) as response:
        html = response.read().decode("utf-8")

    if (html):
        html = html.splitlines()
        fnames = [line.split(' ')[-1] for line in html]

        storm_name = storm_name.upper()

        year = datetime[2:4]
        month = datetime[4:6]
        day = datetime[6:8]
        hour_true = datetime[8:10]

        # Determine the hour for the filename
        if (hour_true not in synop_times):
            if (hour_true < '18'):
                prev_idx = next(idx for idx, value in enumerate(synop_times) if value > hour_true)
                prev_idx -= 1
                hour_adj = synop_times[prev_idx]
            else:
                hour_adj = '18'

            hour_diff = int(hour_true) - int(hour_adj)
            hour = hour_adj
        else:
            hour = hour_true

        ftp_fname = year + month + day + hour + basin + storm_nums[storm_name] + year + '_ships.txt'

        ftp_url = BASE_URL + '/' + ftp_fname

        print("Processing SHIPS file: " + ftp_fname + "\n")

        with urllib.request.urlopen(ftp_url) as response:
            data = response.read().decode("utf-8")

        if (data):
            data = data.splitlines()

            idx = 2
            while (idx < 14):
                line = data[idx].split(' ')
                line = [val for val in line if val != '']

                if (idx == 2):
                    if (storm_name in line):
                        data_dict['datetime'] = datetime
                    else:
                        print("Error: Storm names do not match (ships_parse.fetch_files)")
                        sys.exit(0)
                elif (idx == 10):
                    if (hour_diff == 0):
                        data_dict['shear_spd'] = line[2]
                    else:
                        data_dict['shear_spd'] = line[3]
                elif (idx == 12):
                    if (hour_diff == 0):
                        data_dict['shear_dir'] = line[2]
                    else:
                        data_dict['shear_dir'] = line[3]
                elif (idx == 13):
                    if (hour_diff == 0):
                        data_dict['sst'] = line[2]
                    else:
                        data_dict['sst'] = line[3]
                idx += 1
        else:
            print("Error retrieving SHIPS file data (ships_parse.fetch_files)")
            sys.exit(0)


    else:
        print("Error retrieving SHIPS filenames (ships_parse.fetch_files)")
        sys.exit(0)

    if (write):
        fname = "SHIPS-" + storm_name + "-FullRun.csv"
        dir_path = join(PATH_LINUX_SHIPS, year + "-" + storm_name)

        if not exists(dir_path):
            print("Creating SHIPS subdirectory for " + year + "-" + storm_name + "\n")
            makedirs(dir_path)

        abs_path = join(PATH_LINUX_SHIPS, fname)

        print("Writing dictionary values to csv...\n")
        with open(abs_path, 'a') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(data_dict.values())

    return data_dict



def read_ships_csv(storm_name, year):
    """
    Reads the accumulated SHIPS data csv file

    Parameters
    ----------
    storm_name : str
        Name of the storm being processed
    year : str
        Year the storm being processed occurred


    Returns
    -------
    data_df : Pandas DataFrame
        Pandas dataframe containing accumulated SHIPS data
    """
    col_names = ['datetime', 'shear_spd', 'shear_dir', 'sst']
    dir_path = join(PATH_LINUX_SHIPS, year + "-" + storm_name)

    fname = "SHIPS-" + storm_name + "-FullRun.csv"
    abs_path = join(dir_path, fname)

    if (exists(abs_path)):
        data_df = pd.read_csv(abs_path, sep=",", header=None, dtype=str)
        data_df.columns = col_names

        return data_df

    else:
        print("Error: Local SHIPS file does not exist (ships_parse.read_ships_csv)")
        sys.exit(0)



def df_from_list(data_dict_list, storm_name, write=False):
    """
    Takes a list of data_dict constructed in fetch_files and creates a pandas
    dataframe. Writes the dataframe to csv is 'write' = True

    Parameters
    ----------
    data_dict_list : list of dictionaries of type key: str, value: str
        List of dictionaries holding SHIPS data. Keys: storm, shear_spd,
        shear_dir, sst
    write : bool
        If true, writes the dataframe to a csv. Default: False

    Returns
    -------
    data_df : Pandas dataframe
        Dataframe of accumulated SHIPS data. Column names: datetime, shear_spd,
        shear_dir, sst
    """
    datetime_start = data_dict_list[0]['datetime']
    datetime_end = data_dict_list[-1]['datetime']

    data_df = pd.DataFrame(data_dict_list)

    if (write):
        fname = "SHIPS-" + storm_name + "-" + datetime_start + "-" + datetime_end + ".txt"
        dir_path = join(PATH_LINUX_SHIPS, datetime_start[:4] + "-" + storm_name)

        if not exists(dir_path):
            print("Creating SHIPS subdirectory for " + datetime_start[:4] + "-"
                    + storm_name + "...\n")
            makedirs(dir_path)

        abs_path = join(PATH_LINUX_SHIPS, fname)

        with open(abs_path, 'w') as file:
            data_df.to_csv(file, sep = ',', header=False, index=False)

    return data_df
