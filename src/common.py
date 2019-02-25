"""
Created on Mon Feb 25

@author: Matt Nicholson

Functions common to numerous functions from numerous files in order to avoid
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



def padding_zero(int_to_pad, lim):
    """
    Adds a leading, or padding, zero to an integer if needed. Returns the int
    as a str

    Parameters
    ------------
    int_to_pad : int
        The integer to add padding/leading zero to, if needed

    lim : int
        Used to check to see if the param int needs a padding zero added to it.
        Ex:
        int_to_pad = 11, lim = 10 -> no padding needed, '11' is returned
        int_to_pad = 9, lim = 10 -> padding needed, '09' is returned

    Returns
    ------------
    result : str
        int_to_pad cast as a string, regardless if a padding 0 was added to it
    """
    str_2_pad = str(int_to_pad)
    str_lim = str(lim)
    while (len(str_2_pad) < len(str_lim)):
        str_2_pad = '0' + str_2_pad

    return str_2_pad



def date_time_chunk(start_date_time, end_date_time):
    """
    Splits a start date_time & end date_time into 1-hr chunks

    Ex:
    date_time_chunk('2018091218', '2018091305')

    yields:
    ['201809121800', '201809121900', '201809122000', ... , '201809130400', '201809130500']

    Parameters
    ------------
    start_date_time : str
        Beginning of the date_time chunk. Format: YYYYMMDDHH or YYYYMMDDHHMM

    end_date_time : str
        End of the date_time chunk. Format: YYYYMMDDHH or YYYYMMDDHHMM

    Returns
    ------------
    chunks : list of str
        List of date_times between start_date_time & end_date_time (inclusive),
        split into 1-hr chunks
    """

                    # J   F   M   A   M   J   J   A   S   O   N   D
                    # 1   2   3   4   5   6   7   8   9   10  11  12
                    # 0   1   2   3   4   5   6   7   8   9   10  11
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    leap_years = [2016, 2020, 2024, 2028, 2032]

    chunks = [start_date_time + '00']

    curr_date_time = start_date_time

    while (curr_date_time != end_date_time):
        curr_year = int(curr_date_time[:4])
        curr_month = int(curr_date_time[4:6])
        curr_day = int(curr_date_time[6:8])
        curr_hour = int(curr_date_time[8:10])

        curr_hour += 1
        if (curr_hour > 23):
            curr_hour = 0
            curr_day += 1
            if (curr_day == 31 and days_per_month[curr_month - 1] == 30):
                curr_day = 1
                curr_month += 1
            elif (curr_day == 32 and days_per_month[curr_month - 1] == 31):
                curr_day = 1
                curr_month += 1
            elif (curr_day > 28 and curr_month == 2):
                if (curr_year not in leap_years):
                    curr_day = 1
                    curr_month += 1
                elif (curr_day > 29):
                    curr_day = 1
                    curr_month += 1

        if (curr_month > 12):
            curr_month = 1
            curr_year += 1

        curr_year = str(curr_year)

        curr_month = padding_zero(curr_month, 10)
        curr_day = padding_zero(curr_day, 10)
        curr_hour = padding_zero(curr_hour, 10)

        curr_date_time = curr_year + curr_month + curr_day + curr_hour

        curr_date_time_mins = curr_date_time + '00'

        chunks.append(curr_date_time_mins)

    return chunks
