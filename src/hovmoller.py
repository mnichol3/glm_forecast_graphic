import numpy as np
import matplotlib.pyplot as plt
from glm_tc_graphic import calc_dist
import pandas as pd
import matplotlib
from glm_tc_graphic import quadrant_bounding_box

HIST_PATH_LINUX = '/media/mnichol3/easystore/data/hist/'

def histogram(glm_coords, center_coords):
    dist = []

    glm_lons = glm_coords[0]
    glm_lats = glm_coords[1]

    for idx, lon in enumerate(glm_lons):
        curr_lat = glm_lats[idx]
        curr_dist = calc_dist((lon, curr_lat), center_coords)
        dist.append(curr_dist)

    # Bin edges
    # [  0.  25.  50.  75. 100. 125. 150. 175. 200. 225. 250. 275. 300. 325.
    #    350. 375. 400.]
    hist, bins = np.histogram(dist, bins=16, range=(0,400))

    # [ 12   5   0  11 459 323 452 220 271  82  48  69  44   5  99  27]
    return hist, bins


# TODO: Add bins_fname param to dynamically get bin array length
def hovmoller_plot(hist_fname, quadrant, bins_fname=None):
    """
    Creates a hovmoller plot of lightning flash density

    Parameters
    ----------
    fname : str
        Path & filename of a storm's histogram data file
    """

    data_list = []
    Ys = []
    rmws = []
    Xs = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350,
          375]

    with open(hist_fname, 'r') as f:
        for line in f:
            line = line.rstrip()
            line = line.split(' ')
            Ys.append(line[0][:8] + "-" + line[0][8:] + "z")
            rmws.append(line[1])
            data_list.append(line[2:])

    fig = plt.figure(figsize=(10, 5))
    ax = plt.gca()

    cmap = plt.get_cmap('hot')

    im = plt.contourf(Xs, Ys, data_list, levels=Xs,cmap=cmap)
    ax.yaxis.set_major_locator(plt.MaxNLocator(20))
    plt.title("FLORENCE-2018 " + quadrant + " Quadrant Total Lightning")
    fig.tight_layout()

    plt.show()



def main():
    quadrants = {'RU': 'Right Upshear', 'RD': 'Right Downshear', 'LD': 'Left Downshear',
                 'LU': 'Left Upshear'}

    for key, val in quadrants.items():
        print("Creating Hovmoller plot for Florence " + val + "...\n")
        hist_fname = HIST_PATH_LINUX + 'FLORENCE-2018-' + key +  '.txt'
        data = hovmoller_plot(hist_fname, val)



if __name__ == "__main__":
    main()
