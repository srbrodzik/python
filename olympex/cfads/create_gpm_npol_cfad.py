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
import pdb

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
inDir = '/home/disk/bob/olympex/zebra/moments/npol/rhi/'
outDir = inDir
date = ['20151105','20151112','20151113','20151114','20151115','20151116',
        '20151117','20151118','20151119','20151120','20151121','20151122',
        '20151123','20151124','20151125','20151126','20151130','20151201',
        '20151202','20151203','20151204','20151205','20151206','20151207',
        '20151208','20151209','20151210','20151211','20151212','20151213',
        '20151214','20151215','20151216','20151217','20151218','20151219',
        '20160103','20160104','20160105','20160106','20160108','20160110',
        '20160111','20160112','20160113','20160114','20160115']
#date = ['20151117']

## variable of interest
refl_name = 'DBZ'
refl_mv = -9999.

## cfad inputs
minVal = 0
maxVal = 50
interval = 1
# ------------------------------------- END INPUTS -------------------------------------

# set up logging
log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# determine number of intervals in cfad
numIntervals = (int)((maxVal-minVal)/interval)

firstFile = 1

# process data for each year and month
for i,idate in enumerate(date):
    log.info('i = {} and idate = {}'.format(i,idate) )
    os.chdir(inDir+idate)

    for file in os.listdir(inDir+idate):
        if file.endswith('nc'):

            # determine direction
            if 'east' in file:
                direction = 'east'
            elif 'west' in file:
                direction = 'west'

            log.info('file = {} and direction = {}'.format(file,direction) )

            # open input reflectivity file
            ncid = nc4.Dataset(file, 'r')

            # read vars & if first file, set up cfad arrays
            refl = ncid.variables[refl_name][:]
            if firstFile:
                log.info('   initializing arrays')
                refl_dims = refl.shape
                numLevels = refl_dims[1]
                size2D = refl_dims[2]*refl_dims[3]
                cfadArr_total = np.zeros((numLevels,numIntervals),dtype=np.int)
                cfadArr_east = np.zeros((numLevels,numIntervals),dtype=np.int)
                cfadArr_west = np.zeros((numLevels,numIntervals),dtype=np.int)
                firstFile = 0

            # close file
            ncid.close()
                
            # go through each level
            # foreach non-missing value, increase count in correct bin
            for ilevel in range(0,numLevels):
                refl_level = refl[:,ilevel,:,:]
                refl_valid = refl_level[~refl_level.mask]
                num_valid = refl_valid.shape[0]
                    
                for igrid in range(0,num_valid):
                   
                    refl_value = refl_valid[igrid]
                    #log.info('igrid = {} and refl val = {}'.format(igrid,refl_value) )

                    # 1. make sure all valid refl values are in [minval,maxval] range
                    if refl_value > maxVal:
                        refl_value = maxVal - 1
                    elif refl_value < minVal:
                        refl_value = minVal
                    #log.info('   new refl val = {}'.format(refl_value) )

                    # 2. figure out which interval each refl value is in
                    ival = refl_value - minVal
                    if ival == 0:
                        iint = 0
                    else:
                        iint = (int)(math.floor(ival/interval))

                    # increment cfad arrays
                    cfadArr_total[ilevel,iint] = cfadArr_total[ilevel,iint]+1
                    if direction == 'east':
                        cfadArr_east[ilevel,iint] = cfadArr_east[ilevel,iint]+1
                        #log.info('   increment east cfad' )
                    elif direction == 'west':
                        cfadArr_west[ilevel,iint] = cfadArr_west[ilevel,iint]+1
                        #log.info('   increment west cfad' )

                    #debug
                    #log.info('refl val = {}'.format(refl_value) )
                    #log.info('iint = {}'.format(iint) )
                    #log.info('total sum = {}'.format(np.sum(cfadArr_total)) )
                    #log.info('east  sum = {}'.format(np.sum(cfadArr_east)) )
                    #log.info('west  sum = {}'.format(np.sum(cfadArr_west)) )

            log.info('   total sum = {}'.format(np.sum(cfadArr_total)) )
            log.info('   east  sum = {}'.format(np.sum(cfadArr_east)) )
            log.info('   west  sum = {}'.format(np.sum(cfadArr_west)) )
                        
# output cfad arrays
fid_cfad_total = open(outDir+'cfad_total.txt','w')
fid_cfad_east = open(outDir+'cfad_east.txt','w')
fid_cfad_west = open(outDir+'cfad_west.txt','w')

log.info('BEFORE OUTPUT total sum = {}'.format(np.sum(cfadArr_total)) )
log.info('BEFORE OUTPUT east  sum = {}'.format(np.sum(cfadArr_east)) )
log.info('BEFORE OUTPUT west  sum = {}'.format(np.sum(cfadArr_west)) )

for ilevel in range(0,numLevels):
    
    cfadArr_total[ilevel,:].tofile(fid_cfad_total,sep="\t",format="%d")
    fid_cfad_total.write("\n")
    
    cfadArr_east[ilevel,:].tofile(fid_cfad_east,sep="\t",format="%d")
    fid_cfad_east.write("\n")

    cfadArr_west[ilevel,:].tofile(fid_cfad_west,sep="\t",format="%d")
    fid_cfad_west.write("\n")
    
fid_cfad_total.close()
fid_cfad_east.close()
fid_cfad_west.close()

pdb.set_trace()


# plot cfad arrays

                
