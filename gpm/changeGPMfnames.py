import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt
import os

prefix = 'GPM-2Ku'
suffix = 'nc'
indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/08/ORIG'

for fname in os.listdir(indir):

    if fname.endswith('NAM.nc'):

        print 'fname = ',fname

        # open a netcdf4 file for reading
        ncid = nc4.Dataset(indir+'/'+fname,'r')

        # read base_time variable
        bt = ncid.variables['base_time'][:]

        # close input file
        ncid.close()

        # convert base_time into date and time UTC
        date = dt.datetime.utcfromtimestamp(bt)
        datetime = date.strftime('%Y%m%d.%H%M%S')

        # get info from input fname
        parts = fname.split('_')
        orbit = parts[5]
        fname_out = prefix+'.'+orbit+'.'+datetime+'.'+suffix
        print 'fname_out = ',fname_out

        # rename file
        os.rename(indir+'/'+fname,indir+'/'+fname_out)

