# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:04:02 2019

@author: Matt Nicholson

This file will contain functions that will organize the date times,
gather the glm, abi, & acft data, and make the necessary function calls to
produce the output
"""

import aws_dl


def main():
    #print(aws_dl.get_os() == 'linux')
    #print(aws_dl.date_time_chunk('2018091218', '2018091305'))
    aws_dl.glm_dl(['2018091218', '2018091219'])




if __name__ == "__main__":
    main()
