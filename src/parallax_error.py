# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 14:31:55 2019

@author: Salty Pete
"""

import math



# Returns tan^2(x) based on the trig identity:
#
#     1 - cos(2x)
#    ------------- = tan^2(x)
#     1 + cos(2x)
#
def tan_sqr(deg):
    rad = math.radians(deg)
    num = 1 - math.cos(2 * rad)
    denom = 1 + math.cos(2 * rad)
    result = num / denom
    
    return result



# 4.4
# cos(ups) > 1/P
def calc_viz_lim(P):
    rhs = 1/P
    ups = math.acos(math.radians(rhs))
    return ups
    
    
    
    
# 4.1
# Projection parameter
# Ratio between the distance to sphere center and sphere radius
def calc_proj_param(R, H):
    P = (R + H) / R
    return P



# 4.2
# a = nadir angle
# ups = angular distance from SSP
def calc_nadir_angle(P, ups):
    num = math.sin(math.radians(ups))
    denom = P - math.cos(math.radians(ups))
    tan_a = num / denom
    a = math.atan(math.radians(tan_a))
    return math.radians(a)



# 4.3
def calc_nadir_angle_inv(P, a):
    fp = P - math.sqrt(1 - (P**2 - 1) * tan_sqr(a))
    sp = math.tan(math.radians(a)) / (1 + tan_sqr(a))
    
    sin_ups = fp + sp
    ups = math.asin(math.radians(sin_ups))
    return ups
    


# 4.5
def calc_ssp_dist(lon, lat):
    cos_ups = math.cos(math.radians(lon)) + math.cos(math.radians(lat))
    ups = math.acos(math.radians(cos_ups))
    return math.degrees(ups)



# 4.11
# lil_h = altitude of elevated cloud layers
def corrected_proj_param(big_H, lil_h, R):
    P_h = (R + big_H) / (R + lil_h)
    return P_h
    
    

# 4.12
def lin_parallaxE(ups, lil_h, R, P):
    num = P * math.sin(math.radians(ups))
    denom = P * math.cos(math.radians(ups)) - 1
    delta_ups = (num / denom) * (lil_h / R)
    return delta_ups



# 4.13
# Parallax error as a function of longitude & latitude
def parallaxE_lon_lat(lon, lat, lil_h, P, R):
    rad_lon = math.radians(lon)
    rad_lat = math.radians(lat)
    
    num = math.cos(rad_lon)**2
    num  = num * ((math.cos(rad_lat)**2))
    num = P * math.sqrt(1 - num)
    
    denom = P * (math.cos(math.radians(lon))) * (math.cos(math.radians(lat))) - 1
    
    delta_ups = (num / denom) * (lil_h / R)
    
    return delta_ups



# 4.14
# Longitude & Latitude parallax error linear approximation
def lin_parallaxE_lon_lat(delta_ups, ups, lon, lat):
    delta_lon = lon * (delta_ups / ups)
    delta_lat = lat * (delta_ups / ups)
    
    return (delta_lon, delta_lat)



def parallaxE_dist(delta_ups, R):
    return (delta_ups * R)    

 
    
# Since GOES-16 central longitude is -75.2 deg and not 0
def calc_lon(lon):
    return lon + 75.2



R = 6378    # Earth radius
H = 35786   # Satellite height above Earth surface

cent_lat = 29.93
cent_lon = -71.35
    
P = calc_proj_param(R, H) # = 6.61
ups = calc_ssp_dist(cent_lon, cent_lat)
#delta_ups = lin_parallaxE(30, 12, R, P)
#print(parallaxE_dist(delta_ups,R))
err = parallaxE_lon_lat(calc_lon(-71.35), 29.93, 14, P, R)
print(parallaxE_dist(err, R))
