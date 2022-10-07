# ***** CFAD code for DOW data (masked) *****
# Author: Stacy Brodzik, University of Washington
# Date: December 21, 2016
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

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
inDir = '/home/disk/bob/olympex/zebra/moments/dow_lo_qc_to20km/rhi'
outDir = inDir
date = ['20151106','20151107','20151108','20151109','20151110',
        '20151111','20151112','20151113','20151114','20151115',
        '20151116','20151117','20151118','20151119',
        '20151123','20151124','20151130',
        '20151201','20151202','20151203','20151204','20151205',
        '20151206','20151208','20151209','20151210',
        '20151211','20151212','20151213','20151214',
        '20151217','20151218','20151219',
        '20160104','20160105',
        '20160106',
        '20160111','20160112','20160113','20160115']

## variable of interest
var_name = 'DBZHC_F'
var_mv = -9999.
var_outfile = 'dow_dbzhc_f'

maskFile = '/home/disk/bob/olympex/zeb-data/dow_mask/mask.nc'
mask_var_name = 'mask'

## cfad inputs
min_value = -25
max_value = 50
interval = 0.5

first_cfad_level = 4   # 2km - lowest level = 0km & delta_z = 0.5
# ------------------------------------- END INPUTS -------------------------------------

# set up logging
log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# determine number of intervals in cfad
numIntervals = (int)((max_value-min_value)/interval)

# read mask
ncid = nc4.Dataset(maskFile, 'r')
mask = np.squeeze( ncid.variables[mask_var_name][:] )
ncid.close()

firstFile = True

# process data for each year and month
for i,idate in enumerate(date):
    log.info('i = {} and idate = {}'.format(i,idate) )
    os.chdir(inDir+'/'+idate)

    for file in os.listdir(inDir+'/'+idate):
        if file.endswith('nc'):

            log.info('file = {}'.format(file) )

            # open input reflectivity file
            ncid = nc4.Dataset(file, 'r')

            # read vars & if first file, set up cfad arrays
            var = ncid.variables[var_name][:]
            if firstFile:
                log.info('   initializing arrays')
                var_dims = var.shape
                numLevels = var_dims[1]
                size2D = var_dims[2]*var_dims[3]
                cfadArr_total = np.zeros((numLevels,numIntervals),dtype=np.int)
                firstFile = False

            # close file
            ncid.close()
                
            # go through each level starting at the 2km level (ilevel = 4)
            # foreach non-missing value, increase count in correct bin
            for ilevel in range(first_cfad_level,numLevels):
                var_level = np.squeeze( var[:,ilevel,:,:] )
                choice = np.logical_and( np.greater(mask,0), np.greater(var_level,var_mv) )
                var_valid = np.extract(choice,var_level)
                num_valid = var_valid.shape[0]
                    
                for igrid in range(0,num_valid):
                   
                    var_value = var_valid[igrid]
                    #log.info('igrid = {} and var val = {}'.format(igrid,var_value) )

                    # 1. make sure all valid var values are in [minval,maxval] range
                    #if var_value > max_value:
                    #    var_value = max_value+1
                    #elif var_value < min_value:
                    #    var_value = min_value
                    #log.info('   new var val = {}'.format(var_value) )

                    # 2. figure out which interval each var value is in
                    ival = (int)(math.floor((var_value - min_value)/interval))
                    if ival >= numIntervals:
                        iint = numIntervals - 1
                    elif ival <= 0:
                        iint = 0
                    else:
                        iint = ival

                    if iint == numIntervals-1:
                        log.info('ilevel = {} and var val = {}'.format(ilevel,var_value) )

                    # increment cfad array
                    cfadArr_total[ilevel,iint] = cfadArr_total[ilevel,iint]+1

                    #debug
                    #log.info('var val = {}'.format(var_value) )
                    #log.info('iint = {}'.format(iint) )
                    #log.info('total sum = {}'.format(np.sum(cfadArr_total)) )
                    #log.info('east  sum = {}'.format(np.sum(cfadArr_east)) )
                    #log.info('west  sum = {}'.format(np.sum(cfadArr_west)) )

            log.info('   total sum = {}'.format(np.sum(cfadArr_total)) )
                        
# output cfad array

# AS AN ALTERNATIVE COULD USE np.savetxt(fname,array)

#fid_cfad_total = open(outDir+var_outfile+'_cfad_total.txt','w')
fid_cfad_total = open(outDir+'/'+var_outfile+'_cfad_for_comp_with_npol_to20km.txt','w')

log.info('BEFORE OUTPUT total sum = {}'.format(np.sum(cfadArr_total)) )

for ilevel in range(0,numLevels):
    
    cfadArr_total[ilevel,:].tofile(fid_cfad_total,sep="\t",format="%d")
    fid_cfad_total.write("\n")
        
fid_cfad_total.close()

exit()


# plot cfad arrays

# TESTING STUFF

numLevels = 24
                
fid_cfad_total_norm = open('npol_cfad_total_norm.txt','w')
fid_cfad_east_norm = open('npol_cfad_east_norm.txt','w')
fid_cfad_west_norm = open('npol_cfad_west_norm.txt','w')

for ilevel in range(0,numLevels):
    
    arr_total_norm[ilevel,:].tofile(fid_cfad_total_norm,sep="\t",format="%8.6f")
    fid_cfad_total_norm.write("\n")
    
    arr_east_norm[ilevel,:].tofile(fid_cfad_east_norm,sep="\t",format="%8.6f")
    fid_cfad_east_norm.write("\n")

    arr_west_norm[ilevel,:].tofile(fid_cfad_west_norm,sep="\t",format="%8.6f")
    fid_cfad_west_norm.write("\n")
    
fid_cfad_total_norm.close()
fid_cfad_east_norm.close()
fid_cfad_west_norm.close()
