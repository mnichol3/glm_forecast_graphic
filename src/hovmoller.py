import numpy as np
import matplotlib.pyplot as plt
from glm_tc_graphic import calc_dist

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

    """
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.show()
    """
    # [ 12   5   0  11 459 323 452 220 271  82  48  69  44   5  99  27]
    # hist[0] = 12
    #print(hist)
    return hist
