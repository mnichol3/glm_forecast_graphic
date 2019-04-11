# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:04:02 2019

@author: Matt Nicholson

This file will contain functions that will organize the date times,
gather the glm, abi, & acft data, and make the necessary function calls to
produce the output
"""

#from . import aws_dl
#from . import vortex_data_parse as vdp
from aws_dl import abi_dl, glm_dl
import vortex_data_parse as vdp
import sys
from os import listdir, mkdir
from os.path import isdir, isfile, join
#import utils
from glm_tc_graphic import accumulate_glm_data, plot_mercator, read_file
from common import get_os
import utils
import ships_parse

PATH_LINUX = '/media/mnichol3/easystore/data'
PATH_LINUX_OUT = '/media/mnichol3/easystore/data/imgs'
PATH_WIN = r'D:\Documents\senior-research-data\data'


def get_obs_path(obs_type):
    """
    Determines the proper file path for the observation directory based on the
    operating system & type of observation.

    Parameters
    ----------
    obs_type : str
        Type of observation being processed


    Returns
    -------
    path : str
        Absolute path of the directory that the observation data files are in


    Notes
    -----
    18 Feb 2019
        Only 'vdm' and 'hdob' obs are supported
        Only linux & windows operating systems are supported
    """

    if (obs_type != 'vdm' and obs_type != 'hdob'):
        print('ERROR: Invalid observation type parameter (manager.get_path)')
        sys.exit(0)

    obs_paths = {'vdm' : ['vdm','REPNT2'], 'hdob' : ['hdob', 'AHONT1']}
    os_type = get_os()

    if (os_type == 'linux'):
        base_path = PATH_LINUX
    elif ('win' in os_type):
        base_path = PATH_WIN
    else:
        print('ERROR: Incompatable operating system (manager.get_path)')
        sys.exit(0)

    obs_path = join(base_path, obs_paths[obs_type][0], obs_paths[obs_type][1])

    return obs_path



def get_obs_file(start_date, end_date, storm_name, obs_type, mode):
    """
    Takes the year a storm occured, the name of that storm, and the desired obs
    type, and returns:
        - mode == "a" : A list of all the files in that directory
        - mode == "s" : A list containing the filename of a file that contains
                        both start_date and end_date in the file name string

    Parameters
    ----------
    start_date : str
        Year the storm of interest occured
    start_date : str
        Year the storm of interest occured
    storm_name : str
        Name of the storm of interest
    obs_type : str
        Type of observation file to pull
    mode : str
        String indicating what to return.
        mode == 'a' --> A list containing all the files in the storm observation
                        directory will be returned
        mode == 's' --> A list containing the filename(s) of the file(s) in the
                        storm observation directory that contains both
                        start_date and end_date in the file name string


    Returns
    -------
    files : list of str
        List of strings containing the filenames of the accumulated observation
        data files in the specified storm observation directory

    OR

    f : list of str
        List containing the filename(s) of the file(s) in the storm observation
        directory that contains both start_date and end_date in the file
        name string

    """
    if (type(start_date) != str):
        start_date = str(start_date)
    if (type(end_date) != str):
        end_date = str(end_date)

    year = start_date[:4]

    obs_path = get_obs_path(obs_type)
    subdir = year + '-' + storm_name.upper()
    abs_path = join(obs_path, subdir)

    try:
        files = [(f, abs_path) for f in listdir(abs_path) if isfile(join(abs_path, f))]
    except FileNotFoundError:
        print('Directory not found: ' + subdir)
        print('in manager.get_obs_file')
        sys.exit(0)

    if (mode == 'a'):
        return files
    elif (mode == 's'):
        for f in files:
            if (start_date in f[0] and end_date in f[0]):
                return [f]
    else:
        print('ERROR: Invalid mode parameter (manager.get_obs_file)')
        sys.exit(0)



def get_vdm(start_date, end_date, storm_name):

    try:
        f_info = get_obs_file(start_date, end_date, storm_name, 'vdm', 's')[0]
    except TypeError:
        print("VDM file does not exist locally")
        print("Downloading VDM data files...")
        vdp_df = vdp.vdm_df(start_date, end_date, storm_name)
    else:
        fname = f_info[0]
        fpath = f_info[1]
        f_abs = join(fpath, fname)
        vdp_df = vdp.read_vdm_csv(f_abs)

    return vdp_df



def make_dir(dirs):

    for p in dirs:
        if (not isdir(p)):
            try:
                mkdir(p)
            except OSError:
                print ("Creation of the directory %s failed" % p)
                sys.exit(0)
            else:
                print ("Created the directory %s ..." % p)



def main():
    """
    Notes
    -----

    FLORENCE:
        Meso2 from 201809010900 - 201809101400
        Meso1 from 201809101400 - ?
    """
    # Takes ~35 seconds to produce 1 graphic

    year = '2018'
    storm_name = 'FLORENCE'
    start_date = '201809101400' #'201809082200'
    end_date = '201809140300'

    storm_dict = {'FLORENCE': ['201809010900', '201809140300', 'meso2']}

    subdirs = ['abi', 'glm', 'vdm', 'imgs', 'SHIPS']
    default_octant = "REPNT2"

    print('\nProcessing storm: ' + year + '-' + storm_name + '\n')

    print('Creating data directories...\n')
    for f in subdirs:
        if (f == 'vdm'):
            path1 = join(PATH_LINUX, f, default_octant)
            path2 = join(PATH_LINUX, f, default_octant, year + '-' + storm_name)
            make_dir([path1, path2])
        else:
            path = join(PATH_LINUX, f, year + '-' + storm_name)
            make_dir([path])

    # Get accumulated vdm df
    print('Downloading VDMs...\n')
    vdm_df = get_vdm(start_date, end_date, storm_name)
    coords = vdp.track_interp(vdm_df, year, 'hour')
    #utils.plot_coords_df(coords)

    """
    We have a list of datetimes from coords corresponding to the LPC
    1-hr interpolation.

    Next, we need to accumulate the GLM data
    """

    datetimes = coords['date_time'].tolist()
    datetimes = [n[:-2] for n in datetimes]

    for idx, dt in enumerate(datetimes):

        print('Downloading GLM data for ' + storm_name + '-' + dt + '...\n')
        glm_fnames = glm_dl(dt, storm_name, True)

        print('Downloading ABI data for ' + storm_name + '-' + dt + '...\n')

        if (int(dt) <= 2018091014):
            sector = 'meso2'
        else:
            sector = 'meso1'

        print('ABI sector: ' + sector + '\n')

        abi_fname = abi_dl(dt + '00', sector, band=13)

        print('\nabi fname: ' + abi_fname + '\n')

        print('Filtering GLM data for ' + storm_name + '-' + dt + '...\n')

        curr_row = coords.iloc[idx]

        center_coords = (float(format(curr_row['lons'], '.3f')),
                            float(format(curr_row['lats'], '.3f')))

        rmw = int(curr_row['rmw'])

        glm_data = accumulate_glm_data(dt, center_coords, storm_name)

        print('Parsing ABI data...\n')
        data_dict = read_file(abi_fname)

        print('Retrieving wind shear data...\n')
        ships_data = ships_parse.fetch_file(dt + '00', storm_name, basin='AL', write=True)
        wind_shear = (ships_data['shear_dir'], ships_data['shear_spd'])

        print('Creating graphic for ' + storm_name + '-' + dt + '...\n')
        plot_mercator(data_dict, glm_data, center_coords, rmw, wind_shear, storm_name)

        print('-----------------------------------------------------------------')

        ##### !!! REMOVE !!! #####
        sys.exit(0) # For testing/debugging
        ##########################

if __name__ == "__main__":
    main()
    #print('glm_forecast_graphic: Calling module <manager> as main...')
