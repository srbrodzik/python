# -*- coding: utf-8 -*-
"""
This script calculates the HID and IWP from the gridded CHIVO data
"""

#############
import numpy as np
from netCDF4 import Dataset
from astropy.convolution import convolve, Gaussian2DKernel
import pyart
from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, 
                            csu_dsd, csu_kdp, csu_misc, fundamentals)
                            
import os
from datetime import datetime
import fnmatch


def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    TJL - taken from old Py-ART version
    Calculate Cartesian coordinate from radar coordinates
    Parameters
    ----------
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    az : array
        Azimuth angle of the radar in degrees.
    ele : array
        Elevation angle of the radar in degrees.
    Returns
    -------
    x, y, z : array
        Cartesian coordinates in meters from the radar.
    Notes
    -----
    The calculation for Cartesian coordinate is adapted from equations
    2.28(b) and 2.28(c) of Doviak and Zrnic [1]_ assuming a
    standard atmosphere (4/3 Earth's radius model).
    .. math::
        z = \\sqrt{r^2+R^2+r*R*sin(\\theta_e)} - R
        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})
        x = s * sin(\\theta_a)
        y = s * cos(\\theta_a)
    Where r is the distance from the radar to the center of the gate,
    :math:\\theta_a is the azimuth angle, :math:\\theta_e is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).
    References
    ----------
    .. [1] Doviak and Zrnic, Doppler Radar and Weather Observations, Second
        Edition, 1993, p. 21.
    """
    theta_e = ele * np.pi / 180.0  # elevation angle in radians.
    theta_a = az * np.pi / 180.0  # azimuth angle in radians.
    R = 6371.0 * 1000.0 * 4.0 / 3.0  # effective radius of earth in meters.
    r = rng * 1000.0  # distances to gates in meters.

    z = (r ** 2 + R ** 2 + 2.0 * r * R * np.sin(theta_e)) ** 0.5 - R
    s = R * np.arcsin(r * np.cos(theta_e) / (R + z))  # arc length in m.
    x = s * np.sin(theta_a)
    y = s * np.cos(theta_a)
    return x, y, z


def get_z_from_radar(radar):
    """Input radar object, return z from radar (km, 2D)"""
    azimuth_1D = radar.azimuth['data']
    elevation_1D = radar.elevation['data']
    srange_1D = radar.range['data']
    sr_2d, az_2d = np.meshgrid(srange_1D, azimuth_1D)
    el_2d = np.meshgrid(srange_1D, elevation_1D)[1]
    xx, yy, zz = radar_coords_to_cart(sr_2d/1000.0, az_2d, el_2d)
    return zz + radar.altitude['data']

def check_sounding_for_montonic(sounding):
    """
    So the sounding interpolation doesn't fail, force the sounding to behave
    monotonically so that z always increases. This eliminates data from
    descending balloons.
    """
    snd_T = np.array(sounding.variables['tdry'])
    snd_z = np.array(sounding.variables['alt'])
    
    #snd_T = sounding.soundingdata['temp']  # In old SkewT, was sounding.data
    #snd_z = sounding.soundingdata['hght']  # In old SkewT, was sounding.data
    dummy_z = []
    dummy_T = []
    #if not snd_T.mask[0]: #May cause issue for specific soundings
    dummy_z.append(snd_z[0])
    dummy_T.append(snd_T[0])
    for i, height in enumerate(snd_z):
        if i > 0:
            if snd_z[i] > snd_z[i-1]:# and not snd_T.mask[i]:
                dummy_z.append(snd_z[i])
                dummy_T.append(snd_T[i])
    snd_z = np.array(dummy_z)
    snd_T = np.array(dummy_T)
    return snd_T, snd_z


def interpolate_sounding_to_radar(sounding, radar, zinput):
    """Takes sounding data and interpolates it to every radar gate."""
    #radar_z = radar.variables['z0']
    #radar_z = get_z_from_radar(radar)
    radar_z = zinput*1000
    radar_T = None
    snd_T, snd_z = check_sounding_for_montonic(sounding)
    shape = np.shape(radar_z)
    #rad_z1d = np.array(radar_z)
    rad_z1d = radar_z.ravel()
    rad_T1d = np.interp(rad_z1d, snd_z, snd_T)
    return np.reshape(rad_T1d, shape), radar_z, snd_z, snd_T


def interpolate_sounding_to_radar1(sounding, radar):
    """Takes sounding data and interpolates it to every radar gate."""
    #radar_z = radar.variables['z0']
    radar_z = get_z_from_radar(radar)
    radar_T = None
    snd_T, snd_z = check_sounding_for_montonic(sounding)
    shape = np.shape(radar_z)
    #rad_z1d = np.array(radar_z)
    rad_z1d = radar_z.ravel()
    rad_T1d = np.interp(rad_z1d, snd_z, snd_T)
    return np.reshape(rad_T1d, shape), radar_z

def get_csu_hid(rad_file, snds, files):
    
    data = Dataset(rad_file, 'r')

    x0 = np.array(data.variables['x0'])
    y0 = np.array(data.variables['y0'])
    z0 = np.array(data.variables['z0'])
    dz = np.array(np.squeeze(data.variables['DZ_qc']))
    zdr = np.array(np.squeeze(data.variables['ZDR_qc']))
    kdp = np.array(np.squeeze(data.variables['KDP_qc']))
    rhohv = np.array(np.squeeze(data.variables['RHOHV_qc']))
    y,z,x = np.meshgrid(y0,z0,x0)

    #get the time from the radarfile
    rad_str = rad_file[112:127]
    #rad_str = rad_file[10:25]
    print (rad_str)
    rad_time = datetime.strptime(rad_str, '%Y%m%d_%H%M%S')
    
    #find the nearest time in sounding times
    closest = min(snds, key=lambda sub: abs(sub - rad_time))
    print (closest)
    ind = snds.index(closest)
    
    #extract that file
    snd_file = '/rasmussen-scratch/mrocque/research/relampago/CACTI_ARM_soundings/netcdf/' + files[ind]
    datas = Dataset(snd_file, 'r')
    temp = np.array(datas.variables['tdry'])
    hght = np.array(datas.variables['alt'])

    #get the interpolated temperature and height from sounding
    radar_T, radar_z, snd_z, snd_t = interpolate_sounding_to_radar(datas, data, z)

    #run the HID algorithm with QCed params and temperature from the sounding
    scores = csu_fhc.csu_fhc_summer(dz=dz, zdr=zdr, rho=rhohv, kdp=kdp, use_temp=True, band='C',
                                    T=radar_T, method='linear', use_trap=True)
    scores1 = np.array(scores)
    scores1[dz <= -100] = 0
    
    #get water and ice mass
    mw, mi = csu_liquid_ice_mass.calc_liquid_ice_mass(dz, zdr, z, T=radar_T)
    mw1 = np.array(mw)
    mi1 = np.array(mi)
    
    #sum up ice mass above 5 km to get ice water path
    iwp = np.sum(mi[9:,:,:], axis=0)
    iwp1 = np.array(iwp)

    return (scores1, iwp1, mw1, mi1)


