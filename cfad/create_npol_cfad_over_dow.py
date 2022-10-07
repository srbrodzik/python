# ***** CFAD code for NPOL data *****
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

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
inDir = '/home/disk/bob/olympex/zebra/moments/npol_qc2_ctrOnDow_to20km/rhi'
outDir = inDir
date = ['20151105',
        '20151112','20151113','20151114','20151115',
        '20151116','20151117','20151118','20151119','20151120',
        '20151121','20151122','20151123','20151124','20151125',
        '20151126','20151130',
        '20151201','20151202','20151203','20151204','20151205',
        '20151206','20151207','20151208','20151209','20151210',
        '20151211','20151212','20151214','20151215',
        '20151216','20151217','20151218','20151219',
#        '20160103','20160104','20160105',
        '20160104','20160105',
        '20160106','20160108','20160110',
        '20160111','20160112','20160113','20160114','20160115']

## variable of interest
#var_name = 'DBZ'
var_name = 'DBZ_TOT'
var_mv = -9999.
#var_outfile = 'npol_dbz'
var_outfile = 'npol_dbz_tot'
#var_offset = -4.5
var_offset = -0.0
#var_name = 'ZDR'
#var_mv = -9999.
#var_outfile = 'npol_zdr'

maskFile = '/home/disk/bob/olympex/zeb-data/dow_mask/mask.nc'
mask_var_name = 'mask'

## cfad inputs
min_value = -25
max_value = 50
interval = 0.5
#min_value = 0
#max_value = 50
#interval = 1
#min_value = -5
#max_value = 5
#interval = 0.2

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

        # go through all netcdf files
        if file.endswith('nc'):

            # only look at eastward facing sectors
            if 'east' in file:

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
                    cfadArr_over_dow = np.zeros((numLevels,numIntervals),dtype=np.int)
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

                        # add offset to refl value before binning
                        var_value = var_value + var_offset
                        
                        # 1. make sure all valid var values are in [minval,maxval] range
                        #if var_value > max_value:
                        #    var_value = max_value
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
                        #ival = var_value - min_value
                        #if ival == 0:
                        #    iint = 1
                        #else:
                        #    iint = (int)(math.floor(ival/interval))

                        if iint >= numIntervals-1:
                            log.info('ilevel = {},iint = {}, var val = {}'.format(ilevel,iint,var_value) )

                        # increment cfad array
                        cfadArr_over_dow[ilevel,iint] = cfadArr_over_dow[ilevel,iint]+1
                        #log.info('   increment npol_over_dow cfad' )

                        #debug
                        #log.info('var val = {}'.format(var_value) )
                        #log.info('iint = {}'.format(iint) )
                        #log.info('over_dow  sum = {}'.format(np.sum(cfadArr_over_dow)) )

            log.info('   over_dow  sum = {}'.format(np.sum(cfadArr_over_dow)) )
                        
# output cfad array

# AS AN ALTERNATIVE COULD USE np.savetxt(fname,array)

#fid_cfad_over_dow = open(outDir+'/'+var_outfile+'_cfad_over_dow_to20km.txt','w')
#fid_cfad_over_dow = open(outDir+'/'+var_outfile+'_cfad_over_dow_to20km_adj-4.txt','w')
fid_cfad_over_dow = open(outDir+'/'+var_outfile+'_cfad_over_dow_to20km_adj'+str(var_offset)+'db.txt','w')

log.info('BEFORE OUTPUT over_dow sum = {}'.format(np.sum(cfadArr_over_dow)) )

for ilevel in range(0,numLevels):
    
    cfadArr_over_dow[ilevel,:].tofile(fid_cfad_over_dow,sep="\t",format="%d")
    fid_cfad_over_dow.write("\n")
    
fid_cfad_over_dow.close()

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
