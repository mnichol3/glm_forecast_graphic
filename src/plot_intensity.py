import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from common import get_os, padding_zero, date_time_chunk




def plot_intensity(fname):

    datetimes = []

    df = pd.read_csv(fname, dtype={'Day':str, 'UTC':str, 'N':str,'W':str, 'kt':str})

    date = df['Day'].tolist()
    time = df['UTC'].tolist()

    lat = df['N'].tolist()
    lon = df['W'].tolist()

    lats = [float(x) for x in lat]
    lons = [float(x) * -1 for x in lon]

    date_full = ['201810' + x for x in date]

    wind = df['kt'].tolist()

    wind = [int(x) for x in wind]

    for idx, d in enumerate(date):
        datetimes.append(date_full[idx] + time[idx])


    total_datetimes = date_time_chunk(datetimes[0][:-2], datetimes[-1][:-2])

    total_datetimes = total_datetimes[1:-1]

    datetimes = [int(x) for x in datetimes]
    #total_datetimes = [int(x) for x in total_datetimes]
    """
    idx = 0
    total_datetimes_2 = []
    while (idx < len(total_datetimes)):
        total_datetimes_2.append(total_datetimes[idx])
        total_datetimes_2.append(total_datetimes[idx][:-2] + '30')
        idx += 1
    """
    timespan = np.arange(0, 60, 6)
    timespan_interp = np.arange(0, 53, 0.5)

    wind_interp = np.interp(timespan_interp, timespan, wind)

    #total_datetimes_2 = [int(a) for a in total_datetimes_2]
    #total_datetimes_trimmed = [a[6:8] + '-' + a[8:] + 'z' for a in total_datetimes_2]

    fig, ax = plt.subplots()
    ax.plot(timespan_interp, wind_interp, 'r^-', label = 'Sustained Wind (kt)')
    ax.set_title('Hurricane Michael October 2018 - Max Sustained Winds', fontsize=18)
    ax.set_ylabel('Wind Speed (knots)', color='r', fontsize=16)
    ax.set_xlabel('Date & Time')
    ax.xaxis.set_major_locator(plt.MaxNLocator(20))
    ax.tick_params('y', colors='r')

    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)

    ax.legend(loc = 'upper left', prop={'size': 14})
    plt.tight_layout()
    plt.show()









fname = '/media/mnichol3/easystore/data/intens/MichaelIntensity.csv'
plot_intensity(fname)
