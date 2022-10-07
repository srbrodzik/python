#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:03:13 2019

@author: broneil
"""

import numpy as np

def calculate_dewpoint(T, RH):
    #Tetens Formual from
    # https://andrewsforest.oregonstate.edu/sites/default/files/lter/data/studies/ms01/dewpt_vpd_calculations.pdf

    x = np.log(RH/100) + ((17.269*T) / (237.3 + T))
    Td = (237.3 * x)/(17.269-x)
        
    return Td
#print(calculate_dewpoint(20,20))
