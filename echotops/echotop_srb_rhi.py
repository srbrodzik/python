# Echo top heights

# Assumes input file names are of the form:
# radar.atoll.<radarname>.refl.YYYYMMDD_hhmmss.nc

import numpy as np
import netCDF4 as nc4
import os

#--------------------- START INPUTS ---------------------

inDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rhi_1km_cf'
outDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/echotops_rhi'
nx=300
ny=300
nz=40
dz=0.5

refl_thresh_vals = [0,10,20,30,40,50]
missing_value = -99.

institution = 'University of Washington'
source = 'SPolKa S-band radar'
title = 'Echotop heights at thresholds of 0, 10, 20, 30, 40, 50 dBZ'

## Dates to process
#dates = ['20111001','20111002','20111003','20111004','20111005',
#         '20111006','20111007','20111008','20111009','20111010',
#         '20111011','20111012','20111013','20111014','20111015',
#         '20111016','20111017','20111018','20111019','20111020',
#         '20111021','20111022','20111023','20111024','20111025',
#         '20111026','20111027','20111028','20111029','20111030',
#         '20111031',
#         '20111101','20111102','20111103','20111104','20111105',
#         '20111106','20111107','20111108','20111109','20111110',
#         '20111111','20111112','20111113','20111114','20111115',
#         '20111116','20111117','20111118','20111119','20111120',
#         '20111121','20111122','20111123','20111124','20111125',
#         '20111126','20111127','20111128','20111129','20111130',
#         '20111201','20111202','20111203','20111204','20111205',
#         '20111206','20111207','20111208','20111209','20111210',
#         '20111211','20111212','20111213','20111214','20111215',
#         '20111216','20111217','20111218','20111219','20111220',
#         '20111221','20111222','20111223','20111224','20111225',
#         '20111226','20111227','20111228','20111229','20111230',
#         '20111231',
#         '20120101','20120102','20120103','20120104','20120105',
#         '20120106','20120107','20120108','20120109','20120110',
#         '20120111','20120112','20120113','20120114','20120115',
#         '20120116']
dates = ['20111016']

#----------------------- END INPUTS -----------------------

# number of reflectivity thresholds
nmax = len(refl_thresh_vals)

for date in dates:
    
    print date
    
    for inFile in os.listdir(inDir+'/'+date):
      if inFile.endswith('nc'):
          
        print inFile
  
        #Filename for output 
        outFile = str( outDir+'/'+date+'/'+inFile.replace('refl','echotops') )

        #If output dir does not exist, create it
        dir = os.path.dirname(outFile)
        if not os.path.exists(dir):
          os.makedirs(dir)

        # read in reflectivity data
        ncidin = nc4.Dataset(inDir+'/'+date+'/'+inFile,'r')
        refl = np.array(np.squeeze(ncidin.variables['REFL']))
        ncidin.close()
        echotop = np.zeros((nmax,ny,nx))+missing_value

        # Find echo tops
        for ii in range(nmax):
            maxht = np.zeros((ny,nx))+missing_value
            for nl in range(nz-1,-1,-1):
                maxht = np.where(np.logical_and((refl[nl,:,:] >= refl_thresh_vals[ii]),
                                                (maxht == missing_value)), nl*dz+dz, maxht)
            echotop[ii,:,:] = maxht
            print np.max(maxht)

        # Write echo top values to a file
        ncidout = nc4.Dataset(outFile,'w',format='NETCDF4')
        
        # create dimensions
        refl_thresh_dim = ncidout.createDimension('refl_thresh_dim',nmax) 
        y = ncidout.createDimension('y',ny)
        x = ncidout.createDimension('x',nx)

        # create variables
        refl_thresh = ncidout.createVariable('refl_thresh',np.float32,('refl_thresh_dim'),zlib=True )
        etop = ncidout.createVariable('echo_top',np.float32,('refl_thresh_dim','y','x'),
                                      fill_value=missing_value, zlib=True )

        # create variable attributes
        # refl_thresh
        refl_thresh.long_name = 'reflectivity_threshold'
        refl_thresh.units = 'dBZ'
        # etop
        etop.long_name = 'reflectivity_maximum_height'
        etop.units = 'km'

        # create global variable attributes
        ncidout.institution = institution
        ncidout.source = source
        ncidout.title = title

        # write vars to file
        refl_thresh[:] = refl_thresh_vals
        etop[:,:,:] = echotop

        ncidout.close()
