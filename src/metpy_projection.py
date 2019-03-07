
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import numpy as np
import cartopy.crs as ccrs
import os
import sys
import matplotlib.cm as cm
import xarray as xr


from datetime import datetime, timedelta
from pyproj import Proj
from netCDF4 import Dataset
import metpy

#fname = '/home/mnichol3/Documents/senior-rsch/data/abi/20180912_1457z_Meso1_Ch1.nc'

fname = '/home/mnichol3/Documents/senior-rsch/data/abi/OR_ABI-L2-CMIPM1-M3C13_G16_s20182560119204_e20182560119272_c20182560119314.nc'

target_lat = 31.72
target_lon = -73.41


def get_data(fname):
    data_dict = {}

    f = Dataset(fname)

    scan_start = datetime.strptime(f.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    scan_end = datetime.strptime(f.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')
    file_created = datetime.strptime(f.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')

    # The 't' variable is the scan's midpoint time
    # I'm not a fan of numpy datetime, so I convert it to a regular datetime object

    print('Scan Start    : %s' % scan_start)
    print('Scan End      : %s' % scan_end)
    print('File Created  : %s' % file_created)
    print('Scan Duration : %.2f minutes' % ((scan_end-scan_start).seconds/60))

    #print(C)

    print("%s is %.2f %s" % (f['band_wavelength'].long_name,
                                 f['band_wavelength'][0],
                                 f['band_wavelength'].units))

    data_dict['data'] = f['CMI'][:]

    sat_h = f['goes_imager_projection'].perspective_point_height
    sat_lon = f['goes_imager_projection'].longitude_of_projection_origin
    sat_sweep = f['goes_imager_projection'].sweep_angle_axis

    data_dict['semi_major'] = f['goes_imager_projection'].semi_major_axis
    data_dict['semi_minor'] = f['goes_imager_projection'].semi_minor_axis

    data_dict['scan_start'] = scan_start
    data_dict['scan_end'] = scan_end
    data_dict['file_created'] = file_created
    data_dict['sat_h'] = sat_h
    data_dict['sat_lon'] = sat_lon
    data_dict['sat_sweep'] = sat_sweep


    # The projection x and y coordinates equals the scanning angle (in radians)
    # multiplied by the satellite height See details here:
    # https://proj4.org/operations/projections/geos.html?highlight=geostationary
    data_dict['x'] = f['x'][:] * sat_h
    data_dict['y'] = f['y'][:] * sat_h

    return data_dict



def plot_geos(data_dict):
    x = data_dict['x']
    y = data_dict['y']

    X, Y = np.meshgrid(x,y)

    fig = plt.figure(figsize=(15, 12))

    globe = ccrs.Globe(semimajor_axis=data_dict['semi_major'],
                       semiminor_axis=data_dict['semi_minor'])

    geos = ccrs.Geostationary(central_longitude=data_dict['sat_lon'],
                              satellite_height=data_dict['sat_h'], globe=globe)

    ax = fig.add_subplot(1, 1, 1, projection=geos)

    plt.pcolormesh(X, Y, data_dict['data'], cmap=cm.Greys)

    ax.coastlines(resolution='50m', color='black', linewidth=2)
    ax.add_feature(ccrs.cartopy.feature.STATES)
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(y), max(y))

    ax.scatter(target_lon,target_lat, marker="+", color="r",
               transform=ccrs.PlateCarree(), s = 200, zorder=10)

    plt.title('GOES-16 Test', loc='left', fontweight='semibold', fontsize=15)
    plt.title('%s' % data_dict['scan_start'].strftime('%d %B %Y %H:%M UTC '), loc='right');
    plt.show()



def plot_mercator(data_dict):
    x = data_dict['x']
    y = data_dict['y']
    data = data_dict['data']

    X, Y = np.meshgrid(x,y)

    fig = plt.figure(figsize=(15, 12))

    p = Proj(proj='geos', h=data_dict['sat_h'], lon_0=data_dict['sat_lon'],
             sweep=data_dict['sat_sweep'])

    globe = ccrs.Globe(semimajor_axis=data_dict['semi_major'],
                       semiminor_axis=data_dict['semi_minor'])

    XX, YY = np.meshgrid(x, y)
    lons, lats = p(XX, YY, inverse=True)

    lats[np.isnan(data)] = 57
    lons[np.isnan(data)] = -152

    min_lat = np.min(lats)
    max_lat = np.max(lats)
    min_lon = np.min(lons)

    projection = ccrs.Mercator(central_longitude=data_dict['sat_lon'],
                               min_latitude=min_lat, max_latitude=max_lat,
                               globe=globe, latitude_true_scale=min_lat)

    data_projection = ccrs.PlateCarree()

    ax = plt.axes(projection=projection)

    #ax.set_extent([-71.84, -68.74, 28.48, 30.35], crs=ccrs.PlateCarree())

    newmap = ax.pcolormesh(lons, lats, data, cmap=cm.Greys, linewidth=0, transform=data_projection)

    ax.coastlines(resolution='50m', color='black', linewidth=2)
    ax.add_feature(ccrs.cartopy.feature.STATES)

    ax.scatter(target_lon,target_lat, marker="+", color="r",
               transform=ccrs.PlateCarree(), s = 200, zorder=10)

    plt.title('GOES-16 Test', loc='left', fontweight='semibold', fontsize=15)
    plt.title('%s' % data_dict['scan_start'].strftime('%d %B %Y %H:%M UTC '), loc='right');
    plt.show()




def main():
    data = get_data(fname)
    plot_geos(data)
    #plot_mercator(data)


if __name__ == '__main__':
    main()
