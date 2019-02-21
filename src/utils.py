"""
Various plotting and other utility functions

@author: Matt Nicholson

"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
import numpy as np

BASE_PATH = r'C:\Users\Salty Pete\Desktop\2018_Fall\Senior Research\HurricaneData'



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
