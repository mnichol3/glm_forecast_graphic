# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:06:59 2018

@author: Matt Nicholson
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from .aws_dl import glm_dl, calc_julian_day
#from aws_dl import glm_dl, calc_julian_day
import sys
from os import listdir
from os.path import isfile, join

'''
This script opens & reads downloaded GLM datafiles, accumulating the flash
data into a list containing a list of flash latitudes & a list of
flash longitudes.

file format:
    <sensor>-<level>-<product short name>/<year>/<julian day>/<hour>/OR_<sensor>
    -<level>-<product short name>-M<scanning mode>-C<channel>-G<GOES Satellite>
    -s<start time>_e<end time>_c<central time>.nc

Ex:
data = accumulate_data('2018091214')
plot_data(data)
'''



def get_os():
    """
    Determines the type of operating system being used. Needed for when we are
    loading & saving local files later

    Parameters
    ------------
    none

    Returns
    ------------
    os_type : str
        Type of OS the script is running on
        Ex: 'linux' is the script is running on a Linux OS, such as Ubuntu
    """
    os_type = sys.platform
    return os_type



def accumulate_data(date_time):
    """
    Accumulates GOES-16 GLM data dowloaded from NOAA's Amazon AWS server

    Parameters
    ------------
    date_time : list of str
        Date & time of the desired files, in a 1-hr block. Format: YYYYMMDDHH

    Returns
    ------------
    glm_data : list of str
        List of GLM flash latitudes & longitudes
    """
    flash_lats = np.array([])
    flash_lons = np.array([])
    julian_days = ['a'] * len(date_time)
    times_to_dl = []

    # Make date_time a list if its not one already
    if (type(date_time) == str):
        date_time = [date_time]

    if (get_os() == 'linux'):
        path = '/home/mnichol3/Documents/senior-rsch/data/glm'
    else:
        path = 'D:\Documents\senior-research-data\glm'

    # Get list of already downloaded files
    curr_files = [f for f in listdir(path) if isfile(join(path, f))]

    for x in date_time:
        if (sum(x in fn for fn in curr_files) < 180):
            # Indicates the whole hour has not been downloaded
            times_to_dl.append(x)

    # Only call the download function if needed
    if (times_to_dl != []):
        print('GLM_plotter.py is calling glm_dl!')
        glm_dl(times_to_dl)
    else:
        print('GLM_plotter.py is NOT calling glm_dl!')

    for day in julian_days:
        fnames = [f for f in listdir(path) if isfile(join(path, f)) and day in f]

        for f in fnames:
            file_path = join(path, f)
            fh = Dataset(file_path, mode='r')

            flash_lats = np.append(flash_lats, fh.variables['flash_lat'][:])
            flash_lons = np.append(flash_lons, fh.variables['flash_lon'][:])

            fh.close()

    glm_data = [flash_lons, flash_lats]

    return glm_data



def plot_data(data):
    """
    Creates a simple plot of GLM flash date on a mercator-projection map, with
    each flash represented by a yellow "+"

    Parameters
    ------------
    data : list of lists of str
        List containing 2 list; the first is a list of flash longitude coordinates,
        the second is a list of flash latitude coordinates

    Returns
    ------------
    A plot of GLM flash date on a mercator-projection map, with each flash
    represented by a yellow "+"
    """

    flash_lons = data[0]
    flash_lats = data[1]

    ax = plt.axes(projection=ccrs.Mercator(central_longitude=-70))
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                 name='admin_1_states_provinces_shp')

    ax.add_feature(states, linewidth=.8, edgecolor='gray')

    plt.scatter(flash_lons, flash_lats, marker='+', color='yellow',
                transform=ccrs.PlateCarree())

    plt.show()


if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <GLM_plotter> as main...')
