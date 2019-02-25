"""
Created on Mon Feb 25

@author: Matt Nicholson

Functions common to numerous files in order to avoid
circular import problems
"""

import sys



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
