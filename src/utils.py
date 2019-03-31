"""
Various plotting and other utility functions

@author: Matt Nicholson

"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
import numpy as np
from netCDF4 import Dataset
import boto3
from botocore.handlers import disable_signing
import botocore
from aws_dl import calc_julian_day
import re
import sys

BASE_PATH = r'C:\Users\Salty Pete\Desktop\2018_Fall\Senior Research\HurricaneData'



def read_file():
    """
    Opens & reads a GOES-16 ABI data file, returning a dictionary of data

    !!! NOTE: Returns implroper sat_lon value; return 75.0 but should be 75.2 for
    GOES-16
    """


def plot_coords_tpl(coords):
    """
    Plots geographic coordinates on a map. Used to test observation coordinate
    interpolation fucntions in vortex_data_parse

    Parameters
    ----------
    coords : list of tuples
        List consisting of tuples of (time, lat, lon)
        Coordinates are in decimal degrees

    Returns
    -------
    None - displays a plot of the coordinates
    """

    lats = [x[1] for x in coords]
    lons = [x[2] for x in coords]


    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())

    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='black',
                             name='admin_1_states_provinces_shp')

    land = NaturalEarthFeature('physical', 'land', '50m', facecolor='black')

    ocean = NaturalEarthFeature('physical', 'ocean', '50m', facecolor='black')

    ax.add_feature(land, linewidth=.8, edgecolor='gray', zorder = 1)
    ax.add_feature(states, linewidth=.8, edgecolor='gray', zorder = 2)
    ax.add_feature(ocean, linewidth=.8, edgecolor='gray', zorder = 0)

    ax.plot(lons, lats, transform=ccrs.PlateCarree(), marker = 'o', color = 'r', zorder=3)

    """
    for i, txt in enumerate(dt):
        ax.annotate(txt, (lons[i], lats[i]), color = 'white', zorder = 4)
    """

    plt.axis('equal')
    plt.title('Hurricane Irma Interpolated VDM Low Pressure Center')

    plt.show()



def plot_coords_df(df):
    """
    Plots geographic coordinates on a map.

    Parameters
    ----------
    df : Pandas Dataframe
        Coordinates are in decimal degrees

    Returns
    -------
    None - displays a plot of the coordinates
    """

    lons = df['lons'].tolist()
    lats = df['lats'].tolist()

    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())

    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='black',
                             name='admin_1_states_provinces_shp')

    land = NaturalEarthFeature('physical', 'land', '50m', facecolor='black')

    ocean = NaturalEarthFeature('physical', 'ocean', '50m', facecolor='black')

    ax.add_feature(land, linewidth=.8, edgecolor='gray', zorder = 1)
    ax.add_feature(states, linewidth=.8, edgecolor='gray', zorder = 2)
    ax.add_feature(ocean, linewidth=.8, edgecolor='gray', zorder = 0)

    ax.plot(lons, lats, transform=ccrs.PlateCarree(), marker = 'o', color = 'r', zorder=3)

    """
    for i, txt in enumerate(dt):
        ax.annotate(txt, (lons[i], lats[i]), color = 'white', zorder = 4)
    """

    plt.axis('equal')
    plt.title('Hurricane Irma Interpolated VDM Low Pressure Center')

    plt.show()



def plot_track_lc():
    """
    Plots hurricane track. Gets data from a local file

    Parameters
    ----------

    Returns
    -------
    None - displays a plot of hurricane track
    """

    storms = ['chris', 'florence', 'harvey', 'michael']

    #for x in storms:
    x = storms[2] + '.csv'
    file_path = os.path.join(BASE_PATH, x)

    df = pd.read_csv(file_path, dtype={'date':str, 'lat':str, 'lon':str,
                                       'wind-sstn (kt)':int, 'wind-gst (kt)':int})

    date_time = df['date'].tolist()
    lats = df['lat'].tolist()
    lons = df['lon'].tolist()

    dt = [x[6:] for x in date_time]
    lats = [float(x[:-1]) for x in lats]
    lons = [float(y[:-1])*-1 for y in lons]

    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())

    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='black',
                             name='admin_1_states_provinces_shp')

    land = NaturalEarthFeature('physical', 'land', '50m', facecolor='black')

    ocean = NaturalEarthFeature('physical', 'ocean', '50m', facecolor='black')

    ax.add_feature(land, linewidth=.8, edgecolor='gray', zorder = 1)
    ax.add_feature(states, linewidth=.8, edgecolor='gray', zorder = 2)
    ax.add_feature(ocean, linewidth=.8, edgecolor='gray', zorder = 0)

    ax.plot(lons, lats, transform=ccrs.PlateCarree(), marker = 'o', color = 'r', zorder=3)

    for i, txt in enumerate(dt):
        ax.annotate(txt, (lons[i], lats[i]), color = 'white', zorder = 4)

    plt.axis('equal')
    #plt.title('Hurricane Harvey Track\n0300z 24 Aug 2017 - 2100z 27 Aug 2017')
    plt.show()



def plot_intensity_lc():
    """
    Plots hurricane intensity as a function of wind speed. Gets data from a local
    file

    Parameters
    ----------

    Returns
    -------
    None - displays a plot of hurricane intensity
    """

    storms = ['chris', 'florence', 'harvey', 'michael']

    #for x in storms:
    x = storms[2] + '.csv'
    file_path = os.path.join(BASE_PATH, x)

    df = pd.read_csv(file_path, dtype={'date':str, 'lat':str, 'lon':str,
                                       'wind-sstn (kt)':int, 'wind-gst (kt)':int})

    date_time = df['date'].tolist()
    wind_sustained = df['wind-sstn (kt)'].tolist()
    wind_gust = df['wind-gst (kt)'].tolist()

    dt = [x[6:] for x in date_time]

    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())

    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='black',
                             name='admin_1_states_provinces_shp')

    land = NaturalEarthFeature('physical', 'land', '50m', facecolor='black')

    ocean = NaturalEarthFeature('physical', 'ocean', '50m', facecolor='black')

    ax.add_feature(land, linewidth=.8, edgecolor='gray', zorder = 1)
    ax.add_feature(states, linewidth=.8, edgecolor='gray', zorder = 2)
    ax.add_feature(ocean, linewidth=.8, edgecolor='gray', zorder = 0)

    plt.plot(dt, wind_sustained, label='Sustained')
    plt.plot(dt, wind_gust, label='Gust')

    plt.ylim(0, 160)
    plt.ylabel('Wind (kts)')

    plt.xticks(rotation='vertical')
    plt.xlabel('Date-Time (UTC)')
    plt.grid(True)
    plt.legend()
    plt.show()



def plot_glm_data(data):
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



def explore_netcdf(fname):
    fh = Dataset(fname, mode='r')
    print("Filename: " + fname.split('/')[-1])

    """
    # Display Variables
    for key in fh.variables:
        print(key)
    """

    #print(fh.variables['flash_lon'][0])



    print("Band ID: " + str(fh.variables['band_id'][0]))
    print("Band Wavelength: " + str("%.2f" % fh.variables['band_wavelength'][0]))
    print("Semi-major axis: " + str(fh.variables['goes_imager_projection'].semi_major_axis))
    print("Semi-minor axis: " + str(fh.variables['goes_imager_projection'].semi_minor_axis))
    print("Inverse flattening: " + str(fh.variables['goes_imager_projection'].inverse_flattening))
    print("Latitude of projection origin: " + str(fh.variables['goes_imager_projection'].latitude_of_projection_origin))
    print("Longitude of projection origin: " + str(fh.variables['goes_imager_projection'].longitude_of_projection_origin))





def explore_aws_glm(date):
    """
    Parameters
    ----------
    date : str
        Format : YYYYMMDD
    """
    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    year = date[:4]
    hour = date[-2:]
    julian_day = calc_julian_day(date)
    month = date[4:6]
    day = date[6:8]

    keys = []
    prefix = 'GLM-L2-LCFA/' + year + '/' + julian_day + '/' + hour + '/OR_GLM-L2-LCFA_G16'
    suffix = ''
    kwargs = {'Bucket': 'noaa-goes16', 'Prefix': prefix}

    while True:
        resp = s3.list_objects_v2(**kwargs)
        for obj in resp['Contents']:
            key = obj['Key']
            if key.endswith(suffix):
                keys.append(key)

        try:
            kwargs['ContinuationToken'] = resp['NextContinuationToken']
        except KeyError:
            break

    #print(keys)
    for x in keys:

        fname_match = re.search('e(\d{14})', x) # Matches scan end date time

        if (fname_match):
            fname = fname_match[1] + '.nc'
            local_fname = year + month + day + fname[-8:-3] + '.nc'
            print(fname)
            print(local_fname)



def explore_aws_abi(date_time, sector):
    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    keys = []
    # RadM1 for Meso sector 1
    # RadM2 for Meso sector 2

    if (len(date_time) != 12 ):
        print('ERROR: Invalid date-time string (aws_dl.abi_dl)')
        sys.exit(0)

    if (sector == 'meso1'):
        sector_prefix = 'M1'
    elif (sector == 'meso2'):
        sector_prefix = 'M2'
    elif (sector == 'conus'):
        sector_prefix = 'C'
    else:
        print('Error: Invalid sector parameter!')
        return

    year = date_time[:4]
    month = date_time[4:6]
    day = date_time[6:8]
    hour = date_time[8:10]
    mins = date_time[-2:]
    julian_day = calc_julian_day(date_time)

    prefix = 'ABI-L1b-RadM/' + year + '/' + julian_day + '/' + hour + '/OR_ABI-L1b-Rad' + sector_prefix + '-M3C01_G16'
    suffix = ''
    kwargs = {'Bucket': 'noaa-goes16', 'Prefix': prefix}
    print(prefix)
    while True:
        resp = s3.list_objects_v2(**kwargs)
        for obj in resp['Contents']:
            key = obj['Key']
            if key.endswith(suffix):
                keys.append(key)

        try:
            kwargs['ContinuationToken'] = resp['NextContinuationToken']
        except KeyError:
            break

    print(keys[0])
    """
    for x in keys:

        #fname_match = re.search('e(\d{11})', x) # Matches scan end date time
        fname_match = re.search('e(' + year + julian_day + hour + mins + ')', x)

        if (fname_match):
            print(fname_match[0])
            fname = fname_match[1] + '.nc'
            local_fname = year + month + day + fname[-7:-3] + '.nc'
            print(fname)
            print(local_fname)
    """



def plot_abi(fname):
    fh = Dataset(fname, mode='r')

    radiance = fh.variables['CMI'][:]
    fh.close()
    fh = None

    fig = plt.figure(figsize=(6,6),dpi=200)
    im = plt.imshow(radiance, cmap='Greys')
    cb = fig.colorbar(im, orientation='horizontal')
    cb.set_ticks([1, 100, 200, 300, 400, 500, 600])
    cb.set_label('Radiance (W m-2 sr-1 um-1)')
    plt.show()


def main():
    #file = '/home/mnichol3/Documents/senior-rsch/data/abi/20180912_1457z_Meso1_Ch1.nc'
    #explore_netcdf(file)
    #explore_aws_glm('20180912')
    #explore_aws_abi('201809121257', 'meso1')

    plot_abi('/media/mnichol3/easystore/data/abi/20180908171700.nc')




if __name__ == "__main__":
    main()
