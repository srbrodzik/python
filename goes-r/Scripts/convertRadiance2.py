#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 16:38:13 2019

@author: brodzik
"""
import numpy as np

def radianceToRefl(radiance,kappa):

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
    
    # Apply the formula to convert radiance to reflectance
    reflectance = kappa * radiance
 
    # Make sure all data is in the valid data range
    reflectance = np.maximum(reflectance, 0.0)
    reflectance = np.minimum(reflectance, 1.0)
    
    return reflectance


def radianceToTb(radiance,fk1,fk2,bc1,bc2):

    # Convert Radiance to Brightness Temperature (Ch07-16)

    #See https://www.star.nesdis.noaa.gov/goesr/docs/ATBD/Imagery.pdf p 22 (Table 8)

    #That refererence is out of data.  Use p.28 of this document:
        #https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf
  
    #This conversion can be expressed as: 
    #     Tb = [ fk2/(alog( (fk1/radiance)+ 1)) - bc1 ] / bc2
    #fk1 (planck_fk1), fk2 (planck_fk2), bc1 (planck_bc1) amd bc2 (planck_bc2) are all present in L1b headers
    
    #TO DO: planck constants from file header and use them for conversion to brightness temp
    
    T = ( fk2/(np.log( (fk1/radiance)+ 1)) - bc1 ) / bc2
    
    return T
 
