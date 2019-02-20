# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 14:31:55 2019

@author: Matt Nicholson
"""

import math

R = 6378    # Earth radius
H = 35786   # Satellite height above Earth surface



def tan_sqr(deg):
    """
    Returns tan^2(x), based on the trig identity:

    1 - cos(2x)
    ------------- = tan^2(x)
    1 + cos(2x)

    Parameters
    ------------
    deg : int or float
        Angle, in degrees, to plug into tan^2(x)

    Returns
    ------------
    result : int/float
        tan^2(deg)
    """
    rad = math.radians(deg)
    num = 1 - math.cos(2 * rad)
    denom = 1 + math.cos(2 * rad)
    result = num / denom

    return result



def calc_viz_lim(P):
    """
    Eq. 4.4
    Calculates the visibility limit of the satellite, which is dependent on
    the satellite's altitude.

    Parameters
    ------------
    P : float
        The projection parameter of the satellite

    Returns
    ------------
    ups : float
        The visibility limit of the satellite
    """
    rhs = 1/P
    ups = math.acos(math.radians(rhs))
    return ups



def calc_proj_param(R, H):
    """
    Eq. 4.1
    Calculates the projection parameter of a geostationary satellite, which is
    the ratio between the distance to the earth's center and the
    earth's radius

    Parameters
    ------------
    R : int
        Radius of the the Earth, in km

    H : int
        Satellite's height above the Earth's surface, in km

    Returns
    ------------
    P : float
        Geostationary satellite's projection paramter
    """
    P = (R + H) / R
    return P



def calc_ssp_dist(lon, lat):
    """
    Eq. 4.5
    Calculates the geostationary satellite's Satellite Sub Point (SSP) distance
    ange (ups) for a given longitude & latitude

    Parameters
    ------------
    lon : float
        Longitude coordinate, in decimal degrees

    lat : float
        Latitude coordinate, in decimal degrees

    Returns
    ------------
    ups: float
        The SSP distance angle for the given longitude & latitude coordinates,
        in degrees
    """
    cos_ups = math.cos(math.radians(lon)) + math.cos(math.radians(lat))
    ups = math.acos(math.radians(cos_ups))
    return math.degrees(ups)



def corrected_proj_param(big_H, lil_h, R):
    """
    Eq. 4.11
    Calculates the corrected projection parameter for a cloud target

    Parameters
    ------------
    big_H : int
        Satellite's height above the Earth's surface, in km

    lil_h : int
        Cloud's altitude above the Earth's surface, in km

    R : int
        Earth's radius, in km

    Returns
    ------------
    P_h : float
        Corrected projection parameter
    """
    P_h = (R + big_H) / (R + lil_h)
    return P_h



def lin_parallaxE(ups, lil_h, R, P):
    """
    Eq. 4.12
    Calculates the linearized expression for parallax error as a function of
    the SSP distance angle

    Parameters
    ------------
    ups : float
        The SSP distance angle, in degrees

    lil_h : int/float
        The cloud's altitude above the Earth's surface, in km

    R : int
        Earth's radius, in km

    P : float
        Satellite projection parameter

    Returns
    ------------
    delta_ups : float
        Angular parallex error, in degrees
    """
    num = P * math.sin(math.radians(ups))
    denom = P * math.cos(math.radians(ups)) - 1
    delta_ups = (num / denom) * (lil_h / R)
    return delta_ups



def parallaxE_lon_lat(lon, lat, lil_h, P, R):
    """
    Eq. 13
    Calculates satellite parallax error as a function of the longitude & latitude
    of an elevated target

    Parameters
    ------------
    lon : float
        Longitude coordinate of the elevated target, in decimal degrees

    lat : float
        Latitude coordinate of the elevated target, in decimal degrees

    lil_h : int/float
        Altitude of the elevated target above the Earth's surface

    P : float
        Satellite's projection parameter

    R : int
        Earth's radius, in km

    Returns
    ------------
    delta_ups : float
        Angular parallax error, in degrees
    """
    rad_lon = math.radians(lon)
    rad_lat = math.radians(lat)

    num = math.cos(rad_lon)**2
    num  = num * ((math.cos(rad_lat)**2))
    num = P * math.sqrt(1 - num)

    denom = P * (math.cos(math.radians(lon))) * (math.cos(math.radians(lat))) - 1

    delta_ups = (num / denom) * (lil_h / R)

    return delta_ups



def lin_parallaxE_lon_lat(delta_ups, ups, lon, lat):
    """
    Eq. 4.14
    Linear approximation of the parallax error as a function of the elevated
    target's longitude and latitude

    Parameters
    ------------
    delta_ups : float
        Parallax error, in km

    ups : float
        SSP distance angle, in degrees

    lon : float
        Elevated target's longitude coordinate, in degrees

    lat : float
        Elevated target's latitude coordinate, in degrees

    Returns
    ------------
    tuple containing 2 floats
    [0] = the longitude parallax error
    [1] = the latitude parallax error
    """
    delta_lon = lon * (delta_ups / ups)
    delta_lat = lat * (delta_ups / ups)

    return (delta_lon, delta_lat)



def parallaxE_dist(delta_ups, R):
    """
    Calculates the distance value of the parallax error in km

    Parameters
    ------------
    delta_ups : float
        Angular parallax error, in degrees

    R : int
        Earth's radius, in km

    Returns
    ------------
    delta_ups * R : float
        Parallax error as a distance, in km

    """
    return (delta_ups * R)



def calc_lon(lon):
    """
    Calculates the latitude of an elevated target relative to the GOES-16
    central longitude of -75.2 (75.2 deg W)

    Parameters
    ------------
    lon : float
        Longitude of an elevated target

    Returns
    ------------
    lon + 75.2 : float
        Longitude of the elevated target relative to GOES-16
    """
    return lon + 75.2


"""
cent_lat = 29.93
cent_lon = -71.35

P = calc_proj_param(R, H) # = 6.61
ups = calc_ssp_dist(cent_lon, cent_lat)
#delta_ups = lin_parallaxE(30, 12, R, P)
#print(parallaxE_dist(delta_ups,R))
err = parallaxE_lon_lat(calc_lon(-71.35), 29.93, 14, P, R)
print(parallaxE_dist(err, R))
"""
if __name__ == '__main__':
    print('glm_forecast_graphic: Calling module <parallax_error> as main...')
