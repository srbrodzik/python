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
#from uw_raintype import rtfunctions as func

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
inDir = '/home/disk/bob/gpm/n1_oly_ocean_ku/classify/ex_data/'
outDir = inDir
year = ['2015','2016']
month15 = ['11','12']
month16 = ['01','02','03']

## variable of interest
refl_name = 'refl'
refl_mv = -99.
refl_offset = 1.7
rt_name = 'rain_type'    # 1=stra  2=conv  3=other
rt_mv = -99.

## cfad inputs
#min_value = 0.
#max_value = 50.
min_value = -25.
max_value = 50.
interval = 1

first_cfad_level = 16   # 2km - lowest level = 0km & delta_z = 0.125
# ------------------------------------- END INPUTS -------------------------------------

# set up logging
log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# determine number of intervals in cfad
numIntervals = (int)((max_value-min_value)/interval)

# process data for each year and month
for i,iyear in enumerate(year):
    log.info('iyear = {}'.format(iyear) )
    if iyear == '2015':
        month = month15
    else:
        month = month16
        
    for j,imonth in enumerate(month):
        log.info('imonth = {}'.format(imonth) )
        os.chdir(inDir+iyear+'/'+imonth)
        
        for file in os.listdir(inDir+iyear+'/'+imonth):
            if file.endswith('gz'):
                log.info('file = {}'.format(file) )

                # uncompress file
                os.system('gunzip '+file)
                base = os.path.splitext(file)

                # open input reflectivity file
                ncid = nc4.Dataset(str(base[0]), 'r')

                # read vars & if first file, set up cfad arrays
                refl = ncid.variables[refl_name][:]
                refl_dims = refl.shape
                if i == 0 & j == 0:
                    numLevels = refl_dims[1]
                    cfadArr_total = np.zeros((numLevels,numIntervals),dtype=np.int)
                    cfadArr_stra = np.zeros((numLevels,numIntervals),dtype=np.int)
                    cfadArr_conv = np.zeros((numLevels,numIntervals),dtype=np.int)
                    cfadArr_other = np.zeros((numLevels,numIntervals),dtype=np.int)

                rt = ncid.variables[rt_name][:]
                rt_dims = rt.shape
                size2D = rt_dims[1] * rt_dims[2]

                # close file and recompress
                ncid.close()
                os.system('gzip '+base[0])
                
                # go through easch level starting at the 2km level (ilevel = 16)
                # foreach non-missing value, increase count in correct bin
                for ilevel in range(first_cfad_level,numLevels):
                    refl_level = refl[:,ilevel,:,:]
                    refl_valid = refl_level[~refl_level.mask]
                    rt_valid = rt[~refl_level.mask]
                    num_valid = refl_valid.shape[0]
                    
                    for igrid in range(0,num_valid):
                   
                        refl_value = refl_valid[igrid]
                        
                        rt_value = rt_valid[igrid]

                        # 1. add 1.7dB to refl value before binning
                        #refl_value = func.ZtoDBZ( func.DBZtoZ(refl_value)+1.7 )
                        refl_value = refl_value + refl_offset
                        
                        # 2. make sure all valid refl values are in [minval,maxval] range
                        if refl_value > max_value:
                            refl_value = max_value
                        elif refl_value < min_value:
                            refl_value = min_value
 
                        # 3. figure out which interval each refl value is in
                        ival = refl_value - min_value
                        if ival == 0:
                            iint = 1
                        else:
                            iint = (int)(math.floor(ival/interval))

                        # increment cfad arrays
                        cfadArr_total[ilevel,iint-1] = cfadArr_total[ilevel,iint-1]+1
                        if rt_value == 1:
                            cfadArr_stra[ilevel,iint-1] = cfadArr_stra[ilevel,iint-1]+1
                        elif rt_value == 2:
                            cfadArr_conv[ilevel,iint-1] = cfadArr_conv[ilevel,iint-1]+1
                        elif rt_value == 3:
                            cfadArr_other[ilevel,iint-1] = cfadArr_other[ilevel,iint-1]+1

# output cfad array
fid_cfad_total = open(outDir+'gpm-ku_cfad_total_adj.txt','w')
fid_cfad_stra = open(outDir+'gpm-ku_cfad_stra_adj.txt','w')
fid_cfad_conv = open(outDir+'gpm-ku_cfad_conv_adj.txt','w')
fid_cfad_other = open(outDir+'gpm-ku_cfad_other_adj.txt','w')
for ilevel in range(0,numLevels):
    cfadArr_total[ilevel,:].tofile(fid_cfad_total,sep="\t",format="%d")
    fid_cfad_total.write("\n")
    cfadArr_stra[ilevel,:].tofile(fid_cfad_stra,sep="\t",format="%d")
    fid_cfad_stra.write("\n")
    cfadArr_conv[ilevel,:].tofile(fid_cfad_conv,sep="\t",format="%d")
    fid_cfad_conv.write("\n")
    cfadArr_other[ilevel,:].tofile(fid_cfad_other,sep="\t",format="%d")
    fid_cfad_other.write("\n")
fid_cfad_total.close()
fid_cfad_stra.close()
fid_cfad_conv.close()
fid_cfad_other.close()


# plot cfad arrays

                
