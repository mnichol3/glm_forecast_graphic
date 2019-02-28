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
import aws_dl
import vortex_data_parse as vdp
import sys
from os import listdir
from os.path import isfile, join
#import utils
import glm_tc_graphic
from common import get_os

PATH_LINUX = '/home/mnichol3/Documents/senior-rsch/data/'
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
    files = [(f, abs_path) for f in listdir(abs_path) if isfile(join(abs_path, f))]

    if (mode == 'a'):
        return files
    elif (mode == 's'):
        for f in files:
            if (start_date in f[0] and end_date in f[0]):
                return [f]
    else:
        print('ERROR: Invalid mode parameter (manager.get_obs_file)')
        sys.exit(0)



def main():

    #aws_dl.glm_dl(['2018091218', '2018091218'])
    #glm_data = GLM_plotter.accumulate_data(['2018091218'])
    #GLM_plotter.plot_data(glm_data)


    # vortex_data_parse
    date = '20180913'
    time = '1600'
    storm_name = 'florence'


    #fname = 'VDM-IRMA-201708300000-201709130600.txt'
    #vdm_dict = vdp.get_vdm(date, time, storm_name)
    #print(vdm_dict)
    #vdp.vdm_df('201708300000', '201709130600', 'irma')

    """
    year = '2017'
    storm_name = 'IRMA'
    start_date = '20170830'
    end_date = '20170913'
    f_info = get_obs_file(start_date, end_date, storm_name, 'vdm', 's')[0]
    fname = f_info[0]
    fpath = f_info[1]
    f_abs = join(fpath, fname)
    vdp_df = vdp.read_vdm_csv(f_abs)
    #print(vdp_df['A'].tolist())
    coords = vdp.track_interp(vdp_df, '2018', 'hour')
    #print(coords)
    utils.plot_coords_df(coords)
    """

    """
    data = glm_tc_graphic.accumulate_glm_data('2018091218', (-71.35, 29.93))

    for x in data:
        print(x)
        """

    utils.explore_netcdf('ABI-L2-CMIPM_2018_255_17_OR_ABI-L2-CMIPM1-M3C01_G16_s20182551726204_e20182551726261_c20182551726323.nc')
    #utils.explore_netcdf('/home/mnichol3/Documents/senior-rsch/data/glm/201809121801000.nc')

    #aws_dl.abi_dl('201809121257', 'meso1')
    #print(aws_dl.calc_julian_day('201809121257'))
    #print(vdp.calc_min_list('09032126', '09032245'))
    #print(aws_dl.date_time_chunk('2018091218', '2018091305'))



if __name__ == "__main__":
    main()
    #print('glm_forecast_graphic: Calling module <manager> as main...')
