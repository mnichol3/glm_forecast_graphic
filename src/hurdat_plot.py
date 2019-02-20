# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 16:30:46 2018

@author: Matt Nicholson
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Calculates minutes since 0830 0000z
def mins_since(time_list):
    times = []
    tot_mins = 0
    for x in time_list:

        if (type(x) != str):
            x = str(x)

        curr_time = x[-4:]
        mins = int(curr_time[-2:])
        hrs = int(curr_time[:-2])
        day = int(x[-6:-4])

        # Adjust the day
        if (day >= 30):
            day = day - 30
        else:
            day += 1

        tot_mins = (day * 1440) + (hrs * 60) + mins
        times.append(tot_mins)

    return times







#fname_hurdat = r'C:\Users\Salty Pete\Desktop\2018_Fall\Senior Research\HURDAT2-2017-IRMA.txt'
#fname_vdm = r'C:\Users\Salty Pete\Desktop\2018_Fall\Senior Research\REPNT2\2017-IRMA\VDM-IRMA-201708300000-201709130600.txt'

fname_hurdat = ''
name_vdm =''

def make_plot(fname_hurdat):
    data = pd.read_csv(fname_hurdat, header = None, skiprows = 1, dtype = str)

    # Row 6 = max sustained winds
    # Row 7 = min SLP
    dates = data[0]
    times = data[1]

    dates = [x[4:] for x in dates]

    date_times = []
    i = 0
    while (i < len(dates)):
        curr = dates[i] + times[i][1:]
        date_times.append(curr)
        i += 1

    times = mins_since(date_times)
    max_winds = list(map(int, data[6]))
    min_slp = list(map(int, data[7]))

    vdm = pd.read_csv(fname_vdm, header = None, dtype = str)
    vdm_times = mins_since(vdm[0])
    vdm_mslp = list(map(int, vdm[3]))

    mslp_interp = np.interp(times, vdm_times, vdm_mslp)

    fig, ax1 = plt.subplots(figsize = (13.0, 6.5))
    ax2 = ax1.twinx()
    ax1.plot(times, max_winds, 'r^-', label = 'HURDAT2 Max Sustained Wind')
    ax1.set_title('Hurricane Irma 2017 - HURDAT2 Max Sustained Winds (red) & Min SLP (blue)')
    ax1.set_ylabel('Wind speed (knots)', color='r')
    ax1.set_xlabel('Minutes since 30 Sept 2017 0000z')
    ax1.tick_params('y', colors='r')
    #ax1.set_xticks(times[::2])
    #ax1.set_xticklabels(times[::2], rotation=45, horizontalalignment='right', fontsize=8)

    ax2.plot(times, min_slp, 'bo-', label = 'HURDAT2 Min SLP')
    ax2.plot(times, mslp_interp, 'b^--', label = 'Interp. Acft Obs Min SLP')
    ax2.set_ylabel('Min SLP (mb)', color='b')
    ax2.tick_params('y', colors='b')
    #ax2.set_xticks(times[::2])
    #ax2.set_xticklabels(times[::2], rotation=45, horizontalalignment='right', fontsize=8)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.legend(bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.show()

    plt.savefig('filename1.png')


if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <hurdat_plot> as main...')
