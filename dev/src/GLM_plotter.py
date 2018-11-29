# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:06:59 2018

@author: Salty Pete
"""

from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import pyproj
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as cm
from cartopy.feature import NaturalEarthFeature
from shapely.ops import transform
from shapely.geometry import Point
from matplotlib.patches import Polygon
from functools import partial


'''
file format:
    <sensor>-<level>-<product short name>/<year>/<julian day>/<hour>/OR_<sensor>
    -<level>-<product short name>-M<scanning mode>-C<channel>-G<GOES Satellite>
    -s<start time>_e<end time>_c<central time>.nc
'''

# Calculates the Julian Day date for a given date string. 
# Julian Day date = day since 01 Jan
# 20180912 = 255
def calc_julian_day(date):
    # Format: YYYYMMDDHHMM
    
        
    # Used to calculate julian day
                    # J   F   M   A   M   J   J   A   S   O   N   D
                    # 1   2   3   4   5   6   7   8   9   10  11  12
                    # 0   1   2   3   4   5   6   7   8   9   10  11
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    leap_years = [2016, 2020, 2024, 2028, 2032]
    
    year = date[:4]
    month = date[4:6]
    day = date[6:8]
    
    if (int(year) in leap_years):
        days_per_month[1] += 1
    
    curr_month = int(month) - 2
    julian_day = int(day)
    
    while (curr_month >= 0):
        julian_day += days_per_month[curr_month]
        curr_month -= 1
        
    return julian_day
    


def accumulate_data(fnames):
    event_lats = []
    event_lons = []
    group_lats = []
    group_lons = []
    flash_lats = []
    flash_lons = []
    
    for file in fnames:
        fh = Dataset(file, mode='r')
        
        event_lats.append(fh.variables['event_lat'][:])
        event_lons.append(fh.variables['event_lon'][:])
    
        group_lats.append(fh.variables['group_lat'][:])
        group_lons.append(fh.variables['group_lon'][:])
    
        flash_lats.append(fh.variables['flash_lat'][:])
        flash_lons.append(fh.variables['flash_lon'][:])
    
        fh.close()

    GLM = {'events': [event_lats, event_lons], 'groups': [group_lats, group_lons],
           'flashes': [flash_lats, flash_lons]}
    
    return GLM


def plot_data():
    fname = 'GLM-2018091-1403z.nc'
    file = r"C:\Users\Salty Pete\Desktop\2018_Fall\Senior Research\GLM-2018091-1403z.nc"
    
    fh = Dataset(file, mode='r')
    
    print_flag = False
    
    if (print_flag):
        #print(fh.dimensions.keys())
        print('\nDimensions:\n')
        for x in fh.dimensions.keys():
            print(x)
            
        print('\nVariables:\n')
        for x in fh.variables.keys():
            print(x)
        
        
    event_lat = fh.variables['event_lat'][:]
    event_lon = fh.variables['event_lon'][:]
    
    group_lat = fh.variables['group_lat'][:]
    group_lon = fh.variables['group_lon'][:]
    
    flash_lat = fh.variables['flash_lat'][:]
    flash_lon = fh.variables['flash_lon'][:]
    
    fh.close()
    
    
    ax = plt.axes(projection=ccrs.Mercator(central_longitude=-70))
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                 name='admin_1_states_provinces_shp')
        
    ax.add_feature(states, linewidth=.8, edgecolor='gray')
    
    bins = (800, 895)
    H, y, x = np.histogram2d(event_lon, event_lat, bins=bins)
    Y, X, = np.meshgrid(y,x)
    H = np.transpose(H)
    
    masked = H
    masked = np.ma.array(masked)
    masked[masked == 0] = np.ma.masked
    
    cmesh = plt.pcolormesh(Y, X, masked, transform=ccrs.PlateCarree(), cmap=cm.jet)
    cb = plt.colorbar(orientation='horizontal', pad=.01, shrink=.8, extend='max')
    cb.set_label(r'GLM Flash Count (~14 km$\mathregular{^{2}}$ bins)')
    
    plt.show()
