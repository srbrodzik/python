#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 15:47:37 2019

@author: brodzik
"""
def centerWavelength(channel):

    if channel == '01':
        Center_WL = '(0.47 um)'
    elif channel == '02':
        Center_WL = '(0.64 um)'
    elif channel == '03':
        Center_WL = '(0.865 um)'
    elif channel == '04':
        Center_WL = '(1.378 um)'
    elif channel == '05':
        Center_WL = '(1.61 um)'
    elif channel == '06':
        Center_WL = '(2.25 um)'
    elif channel == '07':
        Center_WL = '(3.90 um)'
    elif channel == '08':
        Center_WL = '(6.19 um)'
    elif channel == '09':
        Center_WL = '(6.95 um)'
    elif channel == '10':
        Center_WL = '(7.34 um)'
    elif channel == '11':
        Center_WL = '(8.50 um)'
    elif channel == '12':
        Center_WL = '(9.61 um)'
    elif channel == '13':
        Center_WL = '(10.35 um)'
    elif channel == '14':
        Center_WL = '(11.20 um)'
    elif channel == '15':
        Center_WL = '(12.30 um)'
    elif channel == '16':
        Center_WL = '(13.30 um)'

    return Center_WL