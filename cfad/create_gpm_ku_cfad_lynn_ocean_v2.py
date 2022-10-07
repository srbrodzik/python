# ***** CFAD code for GPM data *****
# Author: Stacy Brodzik, University of Washington
# Date: September 7, 2016
# Description: 

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import netCDF4 as nc4
import os
import numpy as np
import logging as log
import math
import sys
import gpxpy.geo

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
#inDir = '/home/disk/bob/gpm/south_ca_lynn_ocean_ku/classify/ex_data'
inDir = '/home/disk/bob/gpm/n1_oly_lynn_ocean_ku/classify/ex_data'
outDir = inDir
year = ['2014','2015','2016','2017','2018']
month14 = ['11','12']
month15 = ['01','02','03','11','12']
month16 = ['01','02','03','11','12']
month17 = ['01','02','03','11','12']
month18 = ['01','02','03']

## variable of interest
refl_name = 'refl'
refl_mv = -99.
#refl_offset = 1.7
refl_offset = 0.0
rt_name = 'rain_type'    # 1=stra  2=conv  3=other
rt_mv = -99.

## cfad inputs
#min_value = 0
#max_value = 50
min_value = -25
max_value = 50
interval = 0.5

first_cfad_level = 16

# ------------------------------------- END INPUTS -------------------------------------

# set up logging
log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# determine number of intervals in cfad
numIntervals = (int)((max_value-min_value)/interval)

# process data for each year and month
for i,iyear in enumerate(year):
    log.info('iyear = {}'.format(iyear) )
    #YY = iyear[2:4]
    if iyear == '2014':
        month = month14
    elif iyear == '2015':
        month = month15
    elif iyear == '2016':
        month = month16
    elif iyear == '2017':
        month = month17
    elif iyear == '2018':
        month = month18
    else:
        log.info('iyear = {} not of interest'.format(iyear) )
        
    for j,imonth in enumerate(month):
        log.info('imonth = {}'.format(imonth) )
        os.chdir(inDir+'/'+iyear+'/'+imonth)
        
        for file in os.listdir(inDir+'/'+iyear+'/'+imonth):
            if file.endswith('nc'):
                log.info('file = {}'.format(file) )

                # open input reflectivity file
                ncid = nc4.Dataset(str(file), 'r')

                # read vars & if first file, set up cfad arrays
                refl = ncid.variables[refl_name][:]
                refl_dims = refl.shape
                # remove time dimension - ntime = 1 for all this data
                refl = np.squeeze(refl, axis=(0,))  
                if i == 0 & j == 0:
                    numLevels = refl_dims[1]
                    cfadArr_total = np.zeros((numLevels,numIntervals),dtype=np.int)
                    #cfadArr_stra = np.zeros((numLevels,numIntervals),dtype=np.int)
                    #cfadArr_conv = np.zeros((numLevels,numIntervals),dtype=np.int)
                    #cfadArr_other = np.zeros((numLevels,numIntervals),dtype=np.int)

                rt = ncid.variables[rt_name][:]
                # remove time dimension - ntime = 1 for all this data
                rt = np.squeeze(rt, axis=(0,))
                rt_dims = rt.shape
                size2D = rt_dims[0] * rt_dims[1]

                # close file and recompress
                ncid.close()
                
                # go through easch level
                # foreach non-missing value, increase count in correct bin
                for ilevel in range(first_cfad_level,numLevels):
                    refl_level = refl[ilevel,:,:]
                    refl_valid = refl_level[~refl_level.mask]
                    rt_valid = rt[~refl_level.mask]
                    num_valid = refl_valid.shape[0]
                    
                    for igrid in range(0,num_valid):

                        refl_value = refl_valid[igrid]
                        rt_value = rt_valid[igrid]

                        # add offset to refl value before binning
                        refl_value = refl_value + refl_offset
                        
                        # figure out which bin each refl value is in
                        ival = (int)(math.floor((refl_value - min_value)/interval))
                        if ival >= numIntervals:
                            iint = numIntervals - 1
                        elif ival <= 0:
                            iint = 0
                        else:
                            iint = ival
                        #log.info('iint = {}'.format(iint) )
                        
                        # increment cfad arrays
                        cfadArr_total[ilevel,iint] = cfadArr_total[ilevel,iint]+1
                        #if rt_value == 1:
                        #    cfadArr_stra[ilevel,iint] = cfadArr_stra[ilevel,iint]+1
                        #elif rt_value == 2:
                        #    cfadArr_conv[ilevel,iint] = cfadArr_conv[ilevel,iint]+1
                        #elif rt_value == 3:
                        #    cfadArr_other[ilevel,iint] = cfadArr_other[ilevel,iint]+1

# output cfad arrays
fid_cfad_total = open(outDir+'/gpm-ku_cfad_lynn_ocean_total_NEW2014-18.txt','w')
#fid_cfad_total = open(outDir+'/gpm-ku_cfad_lynn_ocean_total_NEW.txt','w')
#fid_cfad_total = open(outDir+'/gpm-ku_cfad_lynn_ocean_total_adj.txt','w')
#fid_cfad_stra = open(outDir+'/gpm-ku_cfad_lynn_ocean_stra_adj.txt','w')
#fid_cfad_conv = open(outDir+'/gpm-ku_cfad_lynn_ocean_conv_adj.txt','w')
#fid_cfad_other = open(outDir+'/gpm-ku_cfad_lynn_ocean_other_adj.txt','w')
for ilevel in range(0,numLevels):
    cfadArr_total[ilevel,:].tofile(fid_cfad_total,sep="\t",format="%d")
    fid_cfad_total.write("\n")
#    cfadArr_stra[ilevel,:].tofile(fid_cfad_stra,sep="\t",format="%d")
#    fid_cfad_stra.write("\n")
#    cfadArr_conv[ilevel,:].tofile(fid_cfad_conv,sep="\t",format="%d")
#    fid_cfad_conv.write("\n")
#    cfadArr_other[ilevel,:].tofile(fid_cfad_other,sep="\t",format="%d")
#    fid_cfad_other.write("\n")
fid_cfad_total.close()
#fid_cfad_stra.close()
#fid_cfad_conv.close()
#fid_cfad_other.close()

## for testing only
#fid_cfad_total_sub = open(outDir+'/gpm-ku_cfad_lynn_ocean_total_10-50.txt_NEW','w')
#for ilevel in range(0,numLevels):
#    cfadArr_total[ilevel,70:].tofile(fid_cfad_total_sub,sep="\t",format="%d")
#    fid_cfad_total_sub.write("\n")
#fid_cfad_total_sub.close()

# plot cfad arrays

                
