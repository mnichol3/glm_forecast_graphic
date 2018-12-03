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
from aws_bucket import glm_dl 
import os

'''
file format:
    <sensor>-<level>-<product short name>/<year>/<julian day>/<hour>/OR_<sensor>
    -<level>-<product short name>-M<scanning mode>-C<channel>-G<GOES Satellite>
    -s<start time>_e<end time>_c<central time>.nc
'''


def accumulate_data(date_time):
    
    fnames = glm_dl(date_time)
    
    flash_lats = np.array([])
    flash_lons = np.array([])
    
    base_path = 'D:\Documents\senior-research-data\glm'
    
    for file in fnames:
        
        file_path = os.path.join(base_path, file)
        fh = Dataset(file_path, mode='r')
        
        #event_lats.append(fh.variables['event_lat'][:])
        #event_lons.append(fh.variables['event_lon'][:])
    
        #group_lats.append(fh.variables['group_lat'][:])
        #group_lons.append(fh.variables['group_lon'][:])
    
        #flash_lats.append(fh.variables['flash_lat'][:])
        #flash_lons.append(fh.variables['flash_lon'][:])
        
        flash_lats = np.append(flash_lats, fh.variables['flash_lat'][:])
        flash_lons = np.append(flash_lons, fh.variables['flash_lon'][:])
    
        fh.close()

    #GLM = {'events': [event_lats, event_lons], 'groups': [group_lats, group_lons],
           #'flashes': [flash_lats, flash_lons]}
    GLM = [flash_lons, flash_lats]

    return GLM



def plot_data(data):
    
    flash_lons = data[0]
    flash_lats = data[1]
    
    ax = plt.axes(projection=ccrs.Mercator(central_longitude=-70))
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                 name='admin_1_states_provinces_shp')
        
    ax.add_feature(states, linewidth=.8, edgecolor='gray')
    
    plt.scatter(flash_lons, flash_lats, marker='+', color='yellow',
                transform=ccrs.PlateCarree())
    
    plt.show()
    
    
    
# Ex. func calls    
#data = accumulate_data('2018091214')   
#plot_data(data)