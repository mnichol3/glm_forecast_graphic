import numpy as np
import matplotlib.pyplot as plt
import math

# https://stackoverflow.com/questions/37899364/autoscaling-not-working-on-drawing-arrow
def arrow_(ax, plt, x, y, dx, dy, **kwargs):
    ax.arrow(x, y, dx, dy, **kwargs)
    plt.plot([x, x + dx], [y, y + dy], alpha=0)


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
        size=size
    )



def main():
    wind = ('260', '5')
    origin = (2, 3) # origin point

    origin_x = origin[0]
    origin_y = origin[1]

    wind_dir = int(wind[0])
    wind_spd = int(wind[1])

    cartesianAngleRadians = (450-wind_dir)*math.pi/180.0
    terminus_x = origin_x + 5 * math.cos(cartesianAngleRadians) * -1
    terminus_y = origin_y + 5 * math.sin(cartesianAngleRadians) * -1


    ax = plt.axes()
    line = plt.plot([origin_x, terminus_x],[origin_y,terminus_y])[0]
    add_arrow(line, size=30)
    plt.scatter(origin_x, origin_y, marker="+", color="r")
    """
    arrow_(ax, plt, origin[0], origin[1], terminus_x, terminus_y, width=0.05, head_width=0.2,
            head_length=0.1, fc='k', ec='k')
    """
    #ax.arrow(origin[0], origin[1], 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.axis([-10, 10, -10, 10])
    plt.show()



if __name__ == "__main__":
    main()
