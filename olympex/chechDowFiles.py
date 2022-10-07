# Make sure that refl correction fields are present in all DOW data

import netCDF4 as nc4
import numpy as np
import os
import logging as log  

varCorr = ['DBZHCC','DBZHCC_F']

inDir = '/home/disk/bob/olympex/cfradial/moments/dow_lo_qc2/rhi'
outDir = inDir

## Dates to process
#dates = ['20151106','20151107','20151108','20151109','20151110',
#         '20151111','20151112','20151113','20151114','20151115',
#         '20151116','20151117','20151118','20151119',
dates = ['20151119',
         '20151123','20151124',
         '20151130',
         '20151201','20151202','20151203','20151204','20151205',
         '20151206','20151208','20151209','20151210',
         '20151211','20151212','20151213','20151214',
         '20151217','20151218','20151219',
         '20160104','20160105',
         '20160106',
         '20160111','20160112','20160113','20160115']

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for date in dates:
    
  print date
  fileDirDate = inDir+'/'+date
    
  for fname in os.listdir(fileDirDate):
    if fname.endswith('nc'):

      #log.info( "file = {}".format(fname) )
  
      #Open input file
      ncid = nc4.Dataset(str(fileDirDate+'/'+fname),'r')

      #Check to make sure all varCorr are present in input file
      for x in varCorr:
          try:
              ncid.variables[x]
          except:
              log.info( "file = {} missing varCorr".format(fname) )


