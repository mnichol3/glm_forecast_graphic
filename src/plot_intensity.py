# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 18:13:48 2018

@author: Matt Nicholson

Little utility script to plot wind intensity and storm tracks
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature

base_path = r'C:\Users\Salty Pete\Desktop\2018_Fall\Senior Research\HurricaneData'

storms = ['chris', 'florence', 'harvey', 'michael']

#for x in storms:
x = storms[2] + '.csv'
file_path = os.path.join(base_path, x)

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
plt.title('Hurricane Harvey Track\n0300z 24 Aug 2017 - 2100z 27 Aug 2017')

plt.show()
'''
# Plots wind profile
date_time = df['date'].tolist()
wind_sustained = df['wind-sstn (kt)'].tolist()
wind_gust = df['wind-gst (kt)'].tolist()
lats = df['lat'].tolist()
lons = df['lon'].tolist()

dt = [x[6:] for x in date_time]

plt.plot(dt, wind_sustained, label='Sustained')
plt.plot(dt, wind_gust, label='Gust')

plt.ylim(0, 160)
plt.ylabel('Wind (kts)')

plt.xticks(rotation='vertical')
plt.xlabel('Date-Time (UTC)')
plt.grid(True)

plt.title('Hurricane Michael Wind - October 2018')

plt.legend()

plt.show()
'''
