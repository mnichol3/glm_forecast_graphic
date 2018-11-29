# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 14:25:45 2018

@author: Salty Pete
"""



def calc_julian_day(date):
    # Format: YYYYMMDDHHMM
    
        
    # Used to calculate julian day
                    # J   F   M   A   M   J   J   A   S   O   N   D
                    # 1   2   3   4   5   6   7   8   9   10  11  12
                    # 0   1   2   3   4   5   6   7   8   9   10  11
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    leap_years = [2016, 2020, 2024, 2028, 2032]
    
    year = date[:4]
    month = date[4:6]
    day = date[6:8]
    
    if (int(year) in leap_years):
        days_per_month[1] += 1
    
    curr_month = int(month) - 2
    julian_day = int(day)
    
    while (curr_month >= 0):
        julian_day += days_per_month[curr_month]
        curr_month -= 1
        
    print(julian_day)
    

calc_julian_day('20180912')
calc_julian_day('20180914')
