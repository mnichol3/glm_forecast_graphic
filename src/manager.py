# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:04:02 2019

@author: Matt Nicholson

This file will contain functions that will organize the date times,
gather the glm, abi, & acft data, and make the necessary function calls to
produce the output
"""

import aws_dl
import GLM_plotter



def get_os():
    """
    Determines the type of operating system being used. Needed for when we are
    loading & saving local files later

    Parameters
    ------------
    none

    Returns
    ------------
    os_type : str
        Type of OS the script is running on
        Ex: 'linux' is the script is running on a Linux OS, such as Ubuntu
    """
    os_type = sys.platform
    return os_type



def main():
    #print(aws_dl.get_os() == 'linux')
    #print(aws_dl.date_time_chunk('2018091218', '2018091305'))
    #aws_dl.glm_dl(['2018091218', '2018091218'])
    glm_data = GLM_plotter.accumulate_data(['2018091218'])
    #GLM_plotter.plot_data(glm_data)



if __name__ == "__main__":
    main()
