#!/usr/bin/python

# -*- coding: utf-8 -*-
# Paul Hein June 2014
import matplotlib as mp
mp.use('Agg')
import sip
import matplotlib.pyplot as plt
import glob
import argparse
import numpy as np
import netCDF4
import pyart
from copy import deepcopy

# Main Program
if __name__ == '__main__':

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Unfold radar velocities for a day and write the results in CF/Radial format.')
    parser.add_argument('fname', type=str, help='The name of the file')
    #parser.add_argument('startfile', type=str, help='The file to get the start sounding volume from')
    args = parser.parse_args()
    fname = args.fname
    print ('Working on: ' +fname)
    #print ('Using "sounding" file from: ' +args.startfile)
# read in the data
    radar = pyart.io.read(fname,file_field_names=True)
    print radar.fields.keys()
# combine the masks  (noting two Falses is a good point)
    combine = np.ma.mask_or(radar.fields['AZ']['data'].mask,radar.fields['VR']['data'].mask)
    radar.fields['AZ']['data'].mask = np.ma.mask_or(combine, radar.fields['AZ']['data'].mask)
    radar.fields['VR']['data'].mask = np.ma.mask_or(combine,radar.fields['VR']['data'].mask)
    radar.fields['AZ']['data'].data[:]=np.where(combine,radar.fields['AZ']['_FillValue'],radar.fields['AZ']['data'].data)
    radar.fields['VR']['data'].data[:]=np.where(combine,radar.fields['VR']['_FillValue'],radar.fields['VR']['data'].data)
# get startfile (sounding)     
    #sradar = pyart.io.read(args.startfile,file_field_names=True)


# perform dealiasing and plot
    dealias_data=[]
########################## dealias full volume
    #print "Doing 4DD dealias"
    #dealias_data.append(pyart.correct.dealias_fourdd(radar, last_radar=sradar, vel_field = 'VR', corr_vel_field = 'VC', last_vel_field = 'VE', keep_original = True, filt=0))
    #radar.add_field('VC', dealias_data[0])
    #radar.fields['VC']['data'].mask = np.ma.mask_or(combine,radar.fields['VC']['data'].mask)
    print "Doing dealias using region"
    dealias_data.append(pyart.correct.dealias_region_based(radar, refl_field = 'AZ', vel_field = 'VR', corr_vel_field = 'VE', keep_original = True))
########################## delias error
#    print "Full Volume Error in dealias_fourdd has occurred. "

##########################  Full volume save and plot
    print "save data"
    #radar.add_field('VE', dealias_data[2])
    radar.add_field('VE', dealias_data[0])
    radar.fields['VE']['data'].mask = np.ma.mask_or(combine,radar.fields['VE']['data'].mask)
    pyart.io.write_cfradial('start.'+fname[0:-3]+'.unfold.nc',radar,format = 'NETCDF4',time_reference = False)
