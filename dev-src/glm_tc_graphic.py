# -*- coding: utf-8 -*-

# =============================================================================
# Matt Nicholson
# University of Maryland
# Department of Atmospheric & Oceanic Science
# Fall 2018
# 
# This python script recreates the graphical tropical cyclone intensity 
# forecasting aid proposed by Stevenson et al. 2017:
#
# Stevenson, Stephanie N., et al. “A 10-Year Survey of Tropical Cyclone 
# Inner-Core Lightning Bursts and Their Relationship to Intensity Change.” 
# Weather and Forecasting, vol. 33, no. 1, 2018, pp. 23–36., 
# doi:10.1175/waf-d-17-0096.1. 
# 
#
# plot_geos(), georeference(), & plot_mercator() adapted from Brian Blaylock's 
# code:
# https://github.com/blaylockbk/pyBKB_v2/blob/master/BB_goes16/mapping_GOES16_data.ipynb
# 
#
# GOES-16 netCDF files used for testing retrieved from Brian Blaylock's 
# Department of Atmospheric Science, University of Utah website:
# http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
# 
# 
# =============================================================================


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
from GLM_plotter import accumulate_data


###############################################################################
# Opens & reads a GOES-16 ABI data file, returning a dictionary of data
#
# @param    fname (str)                 Name of the GOES-16 ABI data file to be 
#                                       processed
#
# @return   glm_data (list of arrays)   List of GLM flash latitudes & longitudes. 
#                                       
###############################################################################
def read_file(fname):
    
    data_dict = {}
    
    # Ch. 1 - Blue vis. Good resolution, not able to use at night
    #fname = 09142018_2101z_Meso1.nc
    
    # Ch. 11 - Cloud top infrared. Has much lower resolution than visible bands
    #fname = 09142018_2101z_Meso1_ch11.nc
    
    # Ch. 13 - 'Clean' Longwave window. Not much different than Ch. 11
    #fname = 09142018_2101z_Meso1_Ch13.nc
    
    #file = r"C:\Users\Salty Pete\Desktop\2018_Fall\Senior Research\20180912_1457z_Meso1_Ch1.nc"

    fh = Dataset(fname, mode='r')
    
    data_dict['band_id'] = fh.variables['band_id'][0]
    
    if (data_dict['band_id'] < 8):
        print('\n!!! WARNING: Currently plotting non-IR satellite data !!!' )
        
    data_dict['band_wavelength'] = "%.2f" % fh.variables['band_wavelength'][0]
    data_dict['semimajor_ax'] = fh.variables['goes_imager_projection'].semi_major_axis
    data_dict['semiminor_ax'] = fh.variables['goes_imager_projection'].semi_minor_axis
    data_dict['inverse_flattening'] = fh.variables['goes_imager_projection'].inverse_flattening
    data_dict['latitude_of_projection_origin'] = fh.variables['goes_imager_projection'].latitude_of_projection_origin
    data_dict['longitude_of_projection_origin'] = fh.variables['goes_imager_projection'].longitude_of_projection_origin
    data_dict['data_units'] = fh.variables['Rad'].units
    
    # Seconds since 2000-01-01 12:00:00
    add_seconds = fh.variables['t'][0]
    
    # Datetime of scan
    scan_date = datetime(2000, 1, 1, 12) + timedelta(seconds=float(add_seconds))
    
    # Satellite height
    sat_height = fh.variables['goes_imager_projection'].perspective_point_height
    
    # Satellite longitude & latitude
    sat_lon = fh.variables['goes_imager_projection'].longitude_of_projection_origin
    sat_lat = fh.variables['goes_imager_projection'].latitude_of_projection_origin
    
    # Satellite lat/lon extend
    lat_lon_extent = {}
    lat_lon_extent['n'] = fh.variables['geospatial_lat_lon_extent'].geospatial_northbound_latitude
    lat_lon_extent['s'] = fh.variables['geospatial_lat_lon_extent'].geospatial_southbound_latitude
    lat_lon_extent['e'] = fh.variables['geospatial_lat_lon_extent'].geospatial_eastbound_longitude
    lat_lon_extent['w'] = fh.variables['geospatial_lat_lon_extent'].geospatial_westbound_longitude
    
    # Geospatial lat/lon center
    data_dict['lat_center'] = fh.variables['geospatial_lat_lon_extent'].geospatial_lat_center
    data_dict['lon_center'] = fh.variables['geospatial_lat_lon_extent'].geospatial_lon_center
    
    # Satellite sweep
    sat_sweep = fh.variables['goes_imager_projection'].sweep_angle_axis
    
    data = fh.variables['Rad'][:].data
    
    Xs = fh.variables['x'][:]
    Ys = fh.variables['y'][:]
    
    fh.close()
    
    data_dict['scan_date'] = scan_date
    data_dict['sat_height'] = sat_height
    data_dict['sat_lon'] = sat_lon
    data_dict['sat_lat'] = sat_lat
    data_dict['lat_lon_extent'] = lat_lon_extent
    data_dict['sat_sweep'] = sat_sweep
    data_dict['x'] = Xs
    data_dict['y'] = Ys
    data_dict['data'] = data
    
    return data_dict


# =============================================================================
# Plot the GOES-16 netcdf file on a geostationary-projection map
# =============================================================================
def plot_geos(data_dict):

    # The projection x and y coordinates equals
    # the scanning angle (in radians) multiplied by the satellite height 
    # http://proj4.org/projections/geos.html <-- 404'd
    # https://proj4.org/operations/projections/geos.html
    
    sat_height = data_dict['sat_height']
    sat_lon = data_dict['sat_lon']
    sat_sweep = data_dict['sat_sweep']
    scan_date = data_dict['scan_date']
    data = data_dict['data']
    
    Xs = data_dict['x'] * sat_height
    Ys = data_dict['y'] * sat_height
    
    X, Y = np.meshgrid(Xs,Ys)
    
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Geostationary(central_longitude=sat_lon,
                                satellite_height=sat_height,false_easting=0,false_northing=0,
                                globe=None, sweep_axis=sat_sweep))
    
    ax.set_xlim(min(Xs), max(Xs))
    ax.set_ylim(min(Ys), max(Ys))
    
    
    ax.coastlines(resolution='10m', color='gray')
    plt.title('GOES-16 Imagery', fontweight='semibold', fontsize=15)
    plt.title('%s' % scan_date.strftime('%H:%M UTC %d %B %Y'), loc='right')
    plt.pcolormesh(X, Y, data, cmap=cm.Greys_r)
    plt.show()



# =============================================================================
# Calculate latitude & longitude of each point with pyproj
#
# =============================================================================
def georeference(data_dict):
    
    sat_height = data_dict['sat_height']
    sat_lon = data_dict['sat_lon']
    sat_sweep = data_dict['sat_sweep']
    data = data_dict['data'] # (1000, 1000) array
    
    # Multiplying by sat height might not be necessary here
    Xs = data_dict['x'] * sat_height # (1000,)
    Ys = data_dict['y'] * sat_height # (1000,)
    
    p = pyproj.Proj(proj='geos', h=sat_height, lon_0=sat_lon, sweep=sat_sweep)
    
    lons, lats = np.meshgrid(Xs, Ys)
    lons, lats = p(lons, lats, inverse=True)
     
    lats[np.isnan(data)] = np.nan
    lons[np.isnan(data)] = np.nan
    
    return (lons, lats)



# =============================================================================
# Plot the GOES-16 netcdf file on a lambert conformal-projection map
# =============================================================================
def plot_mercator(data_dict):

    # The projection x and y coordinates equals
    # the scanning angle (in radians) multiplied by the satellite height 
    # (http://proj4.org/projections/geos.html)
    
    sat_height = data_dict['sat_height']
    scan_date = data_dict['scan_date']
    data = data_dict['data']
    
    Xs = data_dict['x'] * sat_height
    Ys = data_dict['y'] * sat_height
    
    globe = ccrs.Globe(semimajor_axis=data_dict['semimajor_ax'], semiminor_axis=data_dict['semiminor_ax'], 
                       flattening=None, inverse_flattening=data_dict['inverse_flattening'])
    
    Xs, Ys = georeference(data_dict)
    
    
    fig = plt.figure(figsize=(10, 5))
    
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator(globe=globe))
    
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                             name='admin_1_states_provinces_shp')
    
    ax.add_feature(states, linewidth=.8, edgecolor='black')
    #ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    
    # TODO: For presentation sample, disable title and add it back in on ppt
    plt.title('GOES-16 Ch.' + str(data_dict['band_id']) + ' ' + str(data_dict['band_wavelength']) + 'um', 
              fontweight='semibold', fontsize=12, loc='left')
    print('!!! WARNING: Storm name label currently hard-coded !!!')
    
    plt.title('Florence - %s' % scan_date.strftime('%H:%M UTC %d %B %Y'), loc='right')
    
    # TODO: dynamically get low pressure center lat & lon
    #cent_lon = -78.398732
    #cent_lat = 34.232230
    cent_lat = 29.93
    cent_lon = -71.35
    
    lim_coords = geodesic_point_buffer(cent_lat, cent_lon, 400)
    lats = [float(x[1]) for x in lim_coords.coords[:]]
    lons = [float(x[0]) for x in lim_coords.coords[:]]
    
    min_lon = min(lons)
    max_lon = max(lons)
    
    min_lat = min(lats)
    max_lat = max(lats)
    
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    
#    ax.set_extent([lat_lon_extent['w'], lat_lon_extent['e'], lat_lon_extent['s'], 
#                   lat_lon_extent['n']], crs=ccrs.PlateCarree())
    
    band = data_dict['band_id']
    if (band == 11 or band == 13):
        color = cm.binary
    else:
        color = cm.gray
        
    #color = cm.hsv    
    # cmap hsv looks the coolest
    cmesh = plt.pcolormesh(Xs, Ys, data, transform=ccrs.PlateCarree(), cmap=color)
    
    range_rings = [100, 200, 300]
    
    for x in range_rings:
        
        coord_list = geodesic_point_buffer(cent_lat, cent_lon, x)
        lats = [float(x[1]) for x in coord_list.coords[:]]
        max_lat = max(lats)
       
        # https://stackoverflow.com/questions/27574897/plotting-disconnected-entities-with-shapely-descartes-and-matplotlib
        mpl_poly = Polygon(np.array(coord_list), ec="r", fc="none", transform=ccrs.PlateCarree(), 
                           linewidth=1.25)
        ax.add_patch(mpl_poly)
        plt.text(cent_lon, max_lat + 0.05, str(x) + " km", color = "r", horizontalalignment="center", transform=ccrs.PlateCarree(),
                 fontsize = 15)
    
    # Set lat & lon grid tick marks
    lon_ticks = [x for x in range(-180, 181) if x % 2 == 0]
    lat_ticks = [x for x in range(-90, 91) if x % 2 == 0]
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray', 
                      alpha=0.5, linestyle='--', draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right=False
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    gl.ylabel_style = {'color': 'red', 'weight': 'bold'}
    
    plt.scatter(cent_lon,cent_lat, marker="+", color="r", transform=ccrs.PlateCarree(), 
                s = 200)
    
    glm_data = accumulate_data('2018091214')
    
    flash_lons = glm_data[0]
    flash_lats = glm_data[1]
    
    plt.scatter(flash_lons, flash_lats, marker='+', color='yellow',
                transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(cmesh)
    # Increase font size of colorbar tick labels
    plt.setp(cbar.ax.yaxis.get_ticklabels(), fontsize=12)
    cbar.set_label('Radiance (' + data_dict['data_units'] + ')', fontsize = 14, labelpad = 20)
    plt.tight_layout()
    plt.show()
    


def geodesic_point_buffer(lat, lon, km):
    
    proj_wgs84 = pyproj.Proj(init='epsg:4326')
    
    # Azimuthal equidistant projection
    aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
    project = partial(
        pyproj.transform,
        pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)),
        proj_wgs84)
    buf = Point(0, 0).buffer(km * 1000)  # distance in metres
    
    return transform(project, buf).exterior

    
    
if __name__ == '__main__':
    data_dict = read_file()
    #print(georeference(data_dict))
    plot_mercator(data_dict)
    #plot_geos(data_dict)












