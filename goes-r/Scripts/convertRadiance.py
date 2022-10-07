#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 11:31:41 2019

@author: brodzik
"""

import numpy as np


def radianceToRefl(radiance,channel):

    # Convert Radiance to Reflectance (Ch01-06)

    #Convert Radiance to Reflectance Factor (units range from 0 to 1)
    #It is useful to convert the radiance values from an absolute scale into a 
    #relative scale to make them easier to deal with in the following steps. We 
    #will use the example formula published by NOAA at
    #http://www.goes-r.gov/products/ATBDs/baseline/Imagery_v2.0_no_color.pdf on page 21.

    #That reference is out of date.  Current values for Esun are here:
        #https://www.star.nesdis.noaa.gov/goesr/docs/ATBD/Imagery.pdf on page 21 (Table 6a)
        
    #That refererence is also out of data.  Use pp.27-8 of this document:
        #https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf
  
    #This conversion can be expressed as: reflFactor = kappa * radiance where
    #     kappa = ((pi*d^2)/Esun)
    #kappa (kappa0), d (earth_sun_distance_anomaly_in_AU) and Esun (esun) are all present in L1b headers
    
    #TO DO: read kappa0 from file header and use it for conversion to reflectance factor
    
    # Define some constants needed for the conversion.
    # Esun = in-band solar irradiance at 1 Astronomical Unit (AU)
    # d = the ratio of the actual distance (as acquired via the GRB) to the mean earth-sun distance
    d2 = 0.3
    if channel == '01':
        Esun = 441.868715
    elif channel == '02':
        Esun = 663.274497
    elif channel == '03':
        Esun = 726.721072
    elif channel == '04':
        Esun = 679.143879
    elif channel == '05':
        Esun = 634.454241
    elif channel == '06':
        Esun = 381.148813
    # Apply the formula to convert radiance to reflectance
    reflectance = (radiance * np.pi * d2) / Esun
    # Make sure all data is in the valid data range
    reflectance = np.maximum(reflectance, 0.0)
    reflectance = np.minimum(reflectance, 1.0)
    
    return reflectance


def radianceToTb(radiance,channel):

    # Convert Radiance to Brightness Temperature (Ch07-16)

    #See https://www.star.nesdis.noaa.gov/goesr/docs/ATBD/Imagery.pdf p 22 (Table 8)

    #That refererence is out of data.  Use p.28 of this document:
        #https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf
  
    #This conversion can be expressed as: 
    #     Tb = [ fk2/(alog( (fk1/radiance)+ 1)) - bc1 ] / bc2
    #fk1 (planck_fk1), fk2 (planck_fk2), bc1 (planck_bc1) amd bc2 (planck_bc2) are all present in L1b headers
    
    #TO DO: planck constants from file header and use them for conversion to brightness temp
    
    #Band_07
    if channel == '07':
        fk1 = 2.00774e+05
        fk2 = 3.68909e+03
        bc1 = 0.50777
        bc2 = 0.99929
    elif channel == '08':
        fk1 = 5.03614e+04
        fk2 = 2.32657e+03
        bc1 = 2.12504
        bc2 = 0.99541
    elif channel == '09':
        fk1 = 3.54940e+04
        fk2 = 2.07047e+03 
        bc1 = 0.33291
        bc2 = 0.9992
    elif channel == '10':
        fk1 = 3.00925e+04 
        fk2 = 1.95961e+03
        bc1 = 0.06984
        bc2 = 0.99983
    elif channel == '11':
        fk1 = 1.93733e+04
        fk2 = 1.69207e+03
        bc1 = 0.17462
        bc2 = 0.99951
    elif channel == '12':
        fk1 = 1.34382e+04
        fk2 = 1.49784e+03
        bc1 = 0.10861 
        bc2 = 0.99966
    elif channel == '13':        
        fk1 = 1.07364e+04
        fk2 = 1.38986e+03
        bc1 = 0.13445
        bc2 = 0.99955
    elif channel == '14':
        fk1 = 8.48310e+03
        fk2 = 1.28490e+03
        bc1 = 0.25361
        bc2 = 0.9991
    elif channel == '15':
        fk1 = 6.40146e+03
        fk2 = 1.16980e+03
        bc1 = 0.27049
        bc2 = 0.99894
    elif channel == '16':
        fk1 = 5.06603e+03
        fk2 = 1.08203e+03
        bc1 = 0.07574 
        bc2 = 0.99968
  
    # https://www.star.nesdis.noaa.gov/goesr/docs/ATBD/Imagery.pdf p 22, Equation (3-5)
    T = ( fk2/(np.log( (fk1/radiance)+ 1)) - bc1 ) / bc2
    
    return T
    

        
