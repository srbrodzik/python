#!/usr/bin/python

# Imports and functions
import os
import netCDF4 as nc4
import numpy as np

# User inputs
indir = '/home/disk/bob/gpm/epo_ku/classify/ex_data_v06'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
years = ['2014','2015','2016','2017','2018','2019','2020']

maxFile = ''
maxGridSize = -999
maxNumLats = -999
maxNumLons = -999

for iyear in years:
    for imonth in months:
        if os.path.exists(indir+'/'+iyear+'/'+imonth):
            for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
                if fname.endswith('nc'):

                    print fname

                    # open and read lat and lon from input file
                    ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'r')
                    #ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/outOfBounds/'+fname,'r')
                    nlat = len(ncid.dimensions['lat'])
                    nlon = len(ncid.dimensions['lon'])
                    gridSize = nlat * nlon
                    ncid.close()

                    if gridSize >  maxGridSize:
                        maxFile = fname
                        maxGridSize = gridSize
                        maxNumLats = nlat
                        maxNumLons = nlon

print 'maxFile: ',maxFile,' maxGridSize = ',maxGridSize,' (',maxNumLats,' * ',maxNumLons,')'
