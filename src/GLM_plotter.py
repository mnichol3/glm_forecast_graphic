# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:06:59 2018

@author: Salty Pete
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from aws_dl import glm_dl 
import os

'''
This script opens & reads downloaded GLM datafiles, accumulating the flash
data into a list containing a list of flash latitudes & a list of 
flash longitudes.

file format:
    <sensor>-<level>-<product short name>/<year>/<julian day>/<hour>/OR_<sensor>
    -<level>-<product short name>-M<scanning mode>-C<channel>-G<GOES Satellite>
    -s<start time>_e<end time>_c<central time>.nc
'''

###############################################################################
# Accumulates GOES-16 GLM data downloaded from NOAA's AWS server
#
# @param    date_time (str or list (str))   Date & time of the desired files. Time 
#                                           is 1-hr block. Format: YYYYMMDDHH
#
# @return   glm_data (list of arrays)       List of GLM flash latitudes & longitudes. 
#                                       
###############################################################################
def accumulate_data(date_time):
    
    # Make date_time a list if its not one already
    if (type(date_time) == str):
        date_time = [date_time]
        
    flash_lats = np.array([])
    flash_lons = np.array([])
    
    for x in date_time:
    
        fnames = glm_dl(x)
        
        base_path = 'D:\Documents\senior-research-data\glm'
        
        for file in fnames:
            
            file_path = os.path.join(base_path, file)
            fh = Dataset(file_path, mode='r')
            
            flash_lats = np.append(flash_lats, fh.variables['flash_lat'][:])
            flash_lons = np.append(flash_lons, fh.variables['flash_lon'][:])
        
            fh.close()

    glm_data = [flash_lons, flash_lats]

    return glm_data


# Plots the flash data
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