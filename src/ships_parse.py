# -*- coding: utf-8 -*-
"""
Created on Wed 10 Apr 2019 16:16

@author: Matt Nicholson

These functions obtain and parse Statistical Hurricane Intensity Predictions
Scheme (SHIPS) files

"""
import urllib.request

BASE_URL = 'ftp://ftp.nhc.noaa.gov/atcf/stext'

def fetch_files(date, basin='AL'):
    """
    SHIPS filename format:
        YYMMDDHHBBSSYY_ships.txt
        Where:
            YY : Year
            MM : Month
            DD : Day
            HH : Hour
            BB : Basin (AL = Atlantic, EP = Eastern Pacific)
            SS : Storm number
            YY : Year

        ex:
            18102712AL1618_ships.txt

            27 Oct 2018 12z AL16
    """
    with urllib.request.urlopen(BASE_URL) as response:
        html = response.read().decode("utf-8")
        html = html.splitlines()
        print(html)

        # At this point we have a list of strings containing, among other things,
        # the names of SHIPS files on the FTP server

fetch_files(5)
