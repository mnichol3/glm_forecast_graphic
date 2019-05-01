import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from common import get_os, padding_zero, date_time_chunk
from vortex_data_parse import read_vdm_csv, calc_min_list
import sys



def plot_intensity(fname):
    datetimes = []

    df = pd.read_csv(fname, dtype={'day':str, 'time':str, 'wind':str})

    date = df['day'].tolist()
    time = df['time'].tolist()

    wind = df['wind'].tolist()

    date_full = ['201809' + x for x in date]

    wind = [int(x) for x in wind]

    for idx, d in enumerate(date):
        datetimes.append(date_full[idx] + time[idx])


    total_datetimes = date_time_chunk(datetimes[0][:-2], datetimes[-1][:-2])

    total_datetimes = total_datetimes[1:-1]

    datetimes = [int(x) for x in datetimes]
    #total_datetimes = [int(x) for x in total_datetimes]

    timespan = np.arange(0, 84, 6)
    timespan_interp = np.arange(0, 84, 1)

    wind_interp = np.interp(timespan_interp, timespan, wind)

    print(datetimes)
    sys.exit(0)

    fig, ax = plt.subplots()
    ax.plot(timespan_interp, wind_interp, 'r^-', label = 'Sustained Wind (kt)')
    ax.set_title('Hurricane Florence September 2018 - Max Sustained Winds', fontsize=18)
    ax.set_ylabel('Wind Speed (knots)', color='r', fontsize=16)
    ax.set_xlabel('Date & Time')
    ax.xaxis.set_major_locator(plt.MaxNLocator(20))
    ax.tick_params('y', colors='r')

    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)

    ax.legend(loc = 'upper right', prop={'size': 14})
    plt.tight_layout()
    plt.show()










fname = '/media/mnichol3/easystore/data/intens/FlorenceIntensity.csv'
plot_intensity(fname)
