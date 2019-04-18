# -*- coding: utf-8 -*-

"""
 Matt Nicholson
 University of Maryland
 Department of Atmospheric & Oceanic Science
 Fall 2018

 This python script recreates the graphical tropical cyclone intensity
 forecasting aid proposed by Stevenson et al. 2017:

 Stevenson, Stephanie N., et al. “A 10-Year Survey of Tropical Cyclone
 Inner-Core Lightning Bursts and Their Relationship to Intensity Change.”
 Weather and Forecasting, vol. 33, no. 1, 2018, pp. 23–36.,
 doi:10.1175/waf-d-17-0096.1.

"""


from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import pyproj
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as cm
from cartopy.feature import NaturalEarthFeature
from shapely.ops import transform
from shapely.geometry import Point, LineString
from shapely.geometry import Polygon as shp_Polygon
from shapely.ops import polygonize
from matplotlib.patches import Polygon
from functools import partial
import sys
from os import listdir
from os.path import isfile, join
from math import sin, cos, sqrt, atan2, radians
from common import get_os
import math

PATH_LINUX_ABI = '/media/mnichol3/easystore/data/abi'
PATH_LINUX_GLM = '/media/mnichol3/easystore/data/glm'
PATH_LINUX_OUT = '/media/mnichol3/easystore/data/imgs'

def read_file(abi_file):
    """
    Opens & reads a GOES-16 ABI data file, returning a dictionary of data

    !!! NOTE: Returns implroper sat_lon value; return 75.0 but should be 75.2 for
    GOES-16

    Parameters:
    ------------
    fname : str
        Name of the GOES-16 ABI date file to be opened & processed


    Returns:
    ------------
    data_dict : dictionary of str
        Dictionar of ABI image data & metadata from the netCDF file
    """

    data_dict = {}

    # Ch. 1 - Blue vis. Good resolution, not able to use at night

    # Ch. 11 - Cloud top infrared. Has much lower resolution than visible bands

    # Ch. 13 - 'Clean' Longwave window. Not much different than Ch. 11


    fh = Dataset(join(PATH_LINUX_ABI, abi_file), mode='r')

    data_dict['band_id'] = fh.variables['band_id'][0]

    if (data_dict['band_id'] < 8):
        print('\n!!! WARNING: Currently plotting non-IR satellite data !!!' )

    data_dict['band_wavelength'] = "%.2f" % fh.variables['band_wavelength'][0]
    data_dict['semimajor_ax'] = fh.variables['goes_imager_projection'].semi_major_axis
    data_dict['semiminor_ax'] = fh.variables['goes_imager_projection'].semi_minor_axis
    data_dict['inverse_flattening'] = fh.variables['goes_imager_projection'].inverse_flattening
    data_dict['latitude_of_projection_origin'] = fh.variables['goes_imager_projection'].latitude_of_projection_origin
    data_dict['longitude_of_projection_origin'] = fh.variables['goes_imager_projection'].longitude_of_projection_origin
    data_dict['data_units'] = fh.variables['CMI'].units

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

    data = fh.variables['CMI'][:].data

    Xs = fh.variables['x'][:]
    Ys = fh.variables['y'][:]

    fh.close()
    fh = None

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



def georeference(data_dict):
    """
    Calculates the longitude and latitude coordinates of each data point

    Parameters
    ------------
    data_dict : dictionary
        Dictionary of ABI file data & metadata


    Returns
    ------------
    (lons, lats) : tuple of lists of floats
        Tuple containing a list of data longitude coordinates and a list of
        data latitude coordinates
    """

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

    lats[np.isnan(data)] = 57
    lons[np.isnan(data)] = -152

    return (lons, lats)



def plot_geos(data_dict):
    """
    Plot the GOES-16 ABI file on a geostationary-projection map

    Parameters
    ------------
    data_dict : dictionary
        Dictionary of data & metadata from GOES-16 ABI file


    Returns
    ------------
    A plot of the ABI data on a geostationary-projection map

    The projection x and y coordinates equals the scanning angle (in radians)
    multiplied by the satellite height
    http://proj4.org/projections/geos.html <-- 404'd
    https://proj4.org/operations/projections/geos.html

    """

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

    cent_lat = 29.93
    cent_lon = -71.35

    plt.scatter(cent_lon,cent_lat, marker="+", color="r", transform=ccrs.PlateCarree(),
                s = 200)

    plt.show()



def plot_mercator(data_dict, glm_data, center_coords, rmw, wind_shear, storm_name):
    """
    Plot the GOES-16 data on a lambert-conformal projection map. Includes ABI
    imagery, GLM flash data, 100km, 200km, & 300km range rings, and red "+" at
    the center point

    Parameters
    ------------
    data_dict : dictionary
        Dictionary of data & metadata from GOES-16 ABI file


    Returns
    ------------
    A plot of the ABI data on a geostationary-projection map

    The projection x and y coordinates equals
    the scanning angle (in radians) multiplied by the satellite height
    (http://proj4.org/projections/geos.html)
    """

    scan_date = data_dict['scan_date']
    data = data_dict['data']

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
    plt.title('GOES-16 Ch.' + str(data_dict['band_id']),
              fontweight='semibold', fontsize=10, loc='left')

    plt.title((storm_name + ' - %s\nEnvironmental wind shear: ' + wind_shear[0]
                + '/' + wind_shear[1]) % scan_date.strftime('%H:%M UTC %d %B %Y'),
              fontsize=10, loc='right')

    cent_lat = float(center_coords[1])
    cent_lon = float(center_coords[0])

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

    # Create 400 km range polygon
    coord_list = geodesic_point_buffer(cent_lat, cent_lon, 400)
    bounding_poly = shp_Polygon(coord_list)

    # Uncomment to plot shear-relative quadrant polygons
    # RD = red, RU = green, LU = blue, LD = cyan
    """
    quad_polys = bbox_poly(center_coords, wind_shear[0], bounding_poly)

    poly_colors='rgbc'
    idx = 0
    for key, val in quad_polys.items():
        x,y = val.exterior.coords.xy

        ax.plot(x, y, color=poly_colors[idx], alpha=0.7,
                linewidth=3, solid_capstyle='round', transform=ccrs.PlateCarree())

        idx += 1
    """

    # Draw the RMW
    rmw_ring = geodesic_point_buffer(cent_lat, cent_lon, rmw * 1.852) # convert nm to km
    lats = [float(x[1]) for x in rmw_ring.coords[:]]
    max_lat = max(lats)

    mpl_poly = Polygon(np.array(rmw_ring), ec="fuchsia", fc="none", transform=ccrs.PlateCarree(),
                       linewidth=1.25)
    ax.add_patch(mpl_poly)
    plt.text(cent_lon, max_lat + 0.05, "RMW", color = "fuchsia", horizontalalignment="center", transform=ccrs.PlateCarree(),
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

    # Create shear vector
    shear_vector(ax, plt, cent_lon, cent_lat, wind_shear, width=0.01, head_width=0.08,
            head_length=0.1, fc='k', ec='k', zorder=15)

    plt.scatter(cent_lon,cent_lat, marker="+", color="r", transform=ccrs.PlateCarree(),
                s = 200)

    if (glm_data):
        flash_lons = glm_data[0]
        flash_lats = glm_data[1]

        plt.scatter(flash_lons, flash_lats, marker='+', color='yellow',
                    transform=ccrs.PlateCarree(), zorder=10)
    else:
        print("No GLM flashes within 500 km of current center!\n")

    cbar = plt.colorbar(cmesh,fraction=0.046, pad=0.04)

    # Increase font size of colorbar tick labels
    plt.setp(cbar.ax.yaxis.get_ticklabels(), fontsize=12)
    cbar.set_label('Radiance (' + data_dict['data_units'] + ')', fontsize = 14, labelpad = 20)

    plt.tight_layout()

    fig = plt.gcf()
    fig.set_size_inches((8.5, 11), forward=False)
    fig.savefig(join(PATH_LINUX_OUT, scan_date.strftime('%Y') + '-' + storm_name,
                scan_date.strftime('%Y%m%d-%H%M')) + '.png', dpi=500)

    #plt.show()
    plt.close(fig)

    return bounding_poly



def geodesic_point_buffer(lat, lon, km):
    """
    Creates a circle on on the earth's surface, centered at (lat, lon) with
    radius of km. Used to form the range rings needed for plotting

    Parameters
    ------------
    lat : float
        Latitude coordinate of the circle's center

    lon : float
        Longitude coordinate of the circle's center

    km : int
        Radius of the circle, in km


    Returns
    ------------
    A list of floats that prepresent the coordinates of the circle's edges
    """

    proj_wgs84 = pyproj.Proj(init='epsg:4326')

    # Azimuthal equidistant projection
    aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
    project = partial(
        pyproj.transform,
        pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)),
        proj_wgs84)
    buf = Point(0, 0).buffer(km * 1000)  # distance in metres

    return transform(project, buf).exterior



def calc_dist(coords1, coords2):
    """
    Calculates the distance between a pair of geographic coordinates in decimal-
    degree format

    Parameters
    ----------
    coords1 : Tuple of floats
        coords1[0] = lon 1
        coords1[1] = lat 1

    coords2 : Tuple of floats
        coords2[0] = lon 2
        coords2[1] = lat 2


    Returns
    -------
    dist : float
        Distance between the two coordinates, in km
    """
    R = 6373.0  # Radius of the Earth, in km

    lon1 = radians(coords1[0])
    lat1 = radians(coords1[1])
    lon2 = radians(coords2[0])
    lat2 = radians(coords2[1])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    dist = R * c

    return dist



def accumulate_glm_data(date_time, center_coords, storm_name):
    """
    Accumulates GOES-16 GLM data dowloaded from NOAA's Amazon AWS server and
    saved to a local file

    Parameters
    ------------
    date_time : list of str
        Date & time of the desired files, in a 1-hr block. Format: YYYYMMDDHH

    center_coords : Tuple of float
        Tuple containing the low pressure center's longitude & latitude.
        center_coords[0] = lon
        center_coords[1] = lat

    storm_name : str
        Name of the storm being processed

    Returns
    ------------
    glm_data : list of str
        List of GLM flash latitudes & longitudes

    Notes
    -----


    """
    flash_lats = np.array([])
    flash_lons = np.array([])
    flash_lats_filtered = []
    flash_lons_filtered = []

    if (type(date_time) == str):
        date_time = [date_time]

    #julian_days = ['a'] * len(date_time)
    #times_to_dl = []

    # Make date_time a list if its not one already
    if (type(date_time) == str):
        date_time = [date_time]

    year = date_time[0][:4]

    if (get_os() == 'linux'):
        path = PATH_LINUX_GLM
    else:
        path = 'D:\Documents\senior-research-data\glm'

    path = join(path, year + '-' + storm_name)

    # Get list of already downloaded files
    curr_files = [f for f in listdir(path) if isfile(join(path, f))]

    """
    for x in date_time:
        if (sum(x in fn for fn in curr_files) < 180):
            # Indicates the whole hour has not been downloaded
            times_to_dl.append(x)

    # Only call the download function if needed
    if (times_to_dl != []):
        print('accumulate_glm_data is calling glm_dl...')
        glm_dl(times_to_dl)
    """

    for day in date_time:     # previously : julian_days
        fnames = [f for f in listdir(path) if isfile(join(path, f)) and day in f]
        skipped = 0

        for f in fnames:

            file_path = join(path, f)

            try:
                fh = Dataset(file_path, mode='r')
            except OSError:
                print('!!! WARNING!!! HDF error #' + str(skipped) + ' ' + f)
                if (skipped == 5):
                    print('ERROR: accumulate_glm_data HDF errors exceeded allowable threshold')
                    sys.exit(0)
                else:
                    skipped += 1
            else:
                flash_lats = np.append(flash_lats, fh.variables['flash_lat'][:])
                flash_lons = np.append(flash_lons, fh.variables['flash_lon'][:])

                fh.close()
                fh = None

    #glm_data = list(zip(flash_lons, flash_lats))
    glm_data = (flash_lons, flash_lats)


    # Filter out flashes greater than 500 km away from the low pressure center
    for idx, curr_lon in enumerate(glm_data[0]):
        curr_lat = glm_data[1][idx]

        if (calc_dist((curr_lon, curr_lat), center_coords) < 400.0):
            flash_lons_filtered.append(curr_lon)
            flash_lats_filtered.append(curr_lat)

    return (flash_lons_filtered, flash_lats_filtered)



def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size, transform=ccrs.PlateCarree()
    )




def shear_vector(ax, plt, center_lon, center_lat, wind, **kwargs):
    """
    Calculates and plots the environmental wind shear vector

    Parameters
    ----------
    ax : pyplot axis
        Axis of the pyplot figure the arrow is to be plotted on
    plt : pyplot plot
        Plot of the figure the arrow is to be plotted on
    center_lon : float
        Longitude coordinate of the origin
    center_lat : float
        Latitude coordinate of the origin
    wind : tuple of str
        Tuple containing the wind direction and speed (in kts)
        Format: (direction, speed)

    Returns
    -------
    None (apparently this works)

    """
    wind_dir = int(wind[0])
    wind_spd = int(wind[1])

    cartesianAngleRadians = (450-wind_dir)*math.pi/180.0

    # Multiply by -1 for correct wind vector direction
    terminus_x = center_lon + 2 * math.cos(cartesianAngleRadians) * -1
    terminus_y = center_lat + 2 * math.sin(cartesianAngleRadians) * -1

    line = plt.plot([center_lon, terminus_x],[center_lat,terminus_y], transform=ccrs.PlateCarree())[0]
    #add_arrow(line, size=300)



def quadrant_bounding_box(center_coords, wind_dir):
    angles = [0, 90, 180, 270]
    bbox = []
    r_earth = 6371  # km
    wind_dir = int(wind_dir)

    if (wind_dir <= 180):
        wind_dir += 180
    else:
        wind_dir -= 180

    cent_lon = math.radians(center_coords[0])
    cent_lat = math.radians(center_coords[1])

    #r_dist = math.sqrt(400**2 + 400**2) / r_earth
    r_dist = 450 / r_earth

    for angle in angles:
        curr_bear = wind_dir + angle

        if (curr_bear > 360):
            curr_bear -= 360
        elif (curr_bear < 0):
            curr_bear += 360

        bearing = math.radians(curr_bear)

        """
        r_lat = math.asin(math.sin(cent_lat) * math.cos(r_dist) + math.cos(cent_lat) *
                          math.sin(r_dist) * math.cos(bearing))

        if (math.cos(r_lat) == 0 or abs(math.cos(r_lat)) < 0.000001):
            r_lon = cent_lon
        else:
            r_lon = ((cent_lon - math.asin( math.sin(bearing) * math.sin(r_dist) / math.cos(r_lat))
                     + math.pi ) % (2 * math.pi)) - math.pi
        """

        lat2 = math.asin( math.sin(cent_lat)*math.cos(r_dist) +
                          math.cos(cent_lat)*math.sin(r_dist)*math.cos(bearing))
        lon2 = cent_lon + math.atan2(math.sin(bearing)*math.sin(r_dist)*math.cos(cent_lat),
                                 math.cos(r_dist)-math.sin(cent_lat)*math.sin(lat2))

        lat = math.degrees(lat2)
        lon = math.degrees(lon2)

        bbox.append((lon, lat))

    # [DS right, US right, US left, DS left]
    return bbox



def bbox_poly(center_coords, wind_dir, bounding_poly):
    bbox = quadrant_bounding_box(center_coords, wind_dir)

    parallel_1 = bbox[0]
    parallel_2 = bbox[2]

    perp_1 = bbox[1]
    perp_2 = bbox[3]

    parallel_line = LineString([parallel_1, parallel_2])
    perp_line = LineString([perp_1, perp_2])

    shear_LR = parallel_line.union(LineString(list(bounding_poly.exterior.coords)))

    shear_LR = polygonize(shear_LR)

    shear_L = next(shear_LR)
    shear_R = next(shear_LR)

    """
    x,y = result_R.exterior.coords.xy

    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    ax.plot(x, y, color='#6699cc', alpha=0.7,
            linewidth=3, solid_capstyle='round', zorder=2)
    ax.set_title('Polygon')
    plt.show()
    """
    shear_L_UD = perp_line.union(LineString(list(shear_L.exterior.coords)))
    shear_R_UD = perp_line.union(LineString(list(shear_R.exterior.coords)))

    shear_L_UD = polygonize(shear_L_UD)
    shear_R_UD = polygonize(shear_R_UD)

    shear_L_D = next(shear_L_UD)
    shear_L_U = next(shear_L_UD)

    shear_R_D = next(shear_R_UD)
    shear_R_U = next(shear_R_UD)

    """
    x,y = shear_R_D.exterior.coords.xy

    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    ax.plot(x, y, color='#6699cc', alpha=0.7,
            linewidth=3, solid_capstyle='round', zorder=2)
    ax.set_title('Polygon')
    plt.show()
    """
    quadrant_polys = {'RD': shear_R_D, 'RU': shear_R_U, 'LU': shear_L_U, 'LD': shear_L_D}

    return quadrant_polys



def sort_glm_coords(quadrant_polys, glm_coords):
    """
    Parameters
    ----------
    quadrant_polys : dict of polygons
        Dictionary of polygons representing shear-relative quadrants of the TC
        Keys: RD, RU, LU, LD
    glm_coords : Tuple of lists
        Tuple of list. Format: (flash_lons, flash_lats)

    Returns
    -------
    glm_coords_sorted : Dict of tuples
        Dictionary of tuples containing GLM flash coordinates for a specific
        quadrant. Ex: {'RD': (flash_lons, flash_lats), 'RD': (flash_lons, flash_lats)...}
    """
    glm_coords_sorted = {'RD': ([],[]), 'RU': ([],[]), 'LU': ([],[]), 'LD': ([],[])}
    orphan = 0
    placed = False

    flash_lons = glm_coords[0]
    flash_lats = glm_coords[1]

    for idx, lon in enumerate(flash_lons):
        lat = flash_lats[idx]
        curr_point = Point(lon, lat)

        for quadrant, poly in quadrant_polys.items():
            if (poly.contains(curr_point)):
                glm_coords_sorted[quadrant][0].append(lon)
                glm_coords_sorted[quadrant][1].append(lat)
                placed = True
        if (not placed):
            orphan += 1

    print("Orphan flash coordinates: " + str(orphan) + "\n")

    return glm_coords_sorted



def filter_flash_errors(glm_data, LL_coords, UR_coords):
    """
    There appears to be dense lines of erroneous flashes around 26 N on
    10 Sep & 11 Sep from 1500-2100z. This function will remove these, but is unable
    to distinguish if a flash is genuine or erroneous.

    Parameters
    ----------
    glm_data : list of str
        List of GLM flash latitudes & longitudes
    LL_coords : tuple of str
        Lower lefthand coordinates of the bounding box contaning the area
        of false flashes
    UR_coords : tuple of str
        Upper righthand coordinates of the bounding box contaning the area
        of false flashes

    Returns
    -------
    filtered_flashes : tuple of lists
        Filtered GLM flash coordinates. Format: (flash_lons, flash_lats)
    """
    filtered_lons = []
    filtered_lats = []

    lons = glm_data[0]
    lats = glm_data[1]

    min_lat = LL_coords[1]
    max_lat = UR_coords[1]

    min_lon = LL_coords[0]
    max_lon = UR_coords[0]

    for idx, lon in enumerate(lons):
        lat = lats[idx]

        if ((lat < min_lat or lat > max_lat) or (lon < min_lon or lon > max_lon)):
            filtered_lons.append(lon)
            filtered_lats.append(lat)

    return (filtered_lons, filtered_lats)




if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <glm_tc_graphic> as main...')
