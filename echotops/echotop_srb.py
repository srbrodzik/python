# Echo top heights

# Assumes input file names are of the form:
# radar.atoll.<radarname>.refl.YYYYMMDD_hhmmss.nc

import numpy as np
import netCDF4 as nc4
import time as tm
import datetime as dt
import os

#--------------------- START INPUTS ---------------------

inDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf/refl'
outDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/echotops_sur'
#nx=300
#ny=300
#nz=40
#dz=0.5

refl_thresh_vals = [0,10,20,30,40,50]
missing_value = -99.

Conventions = 'CF-1.0'
title = 'Echotop heights at thresholds of 0, 10, 20, 30, 40, 50 dBZ'
institution = 'University of Washington'
source = 'SPolKa radar data'
references = ''
comment = 'Based on interpolated sur vols created from polar vols with elev angles at 0.5,1.5,2.5,3.5,5,7,9,11 degrees'

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
dates = ['20111025','20111026']

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

        # open radar file
        ncidin = nc4.Dataset(inDir+'/'+date+'/'+inFile,'r')

        # read in dimensions
        nx = len(ncidin.dimensions['x0'])
        ny = len(ncidin.dimensions['y0'])
        nz = len(ncidin.dimensions['z0'])
        nt = len(ncidin.dimensions['time'])
          
        # read in reflectivity data
        refl = np.array(np.squeeze(ncidin.variables['REFL']))

        # read in ancilliary data and attributes
        time_val = np.array(np.squeeze(ncidin.variables['time']))
        date1 = dt.datetime.utcfromtimestamp(time_val)
        datetime = date1.strftime('%Y-%m-%dT%H:%M:%SZ')
        
        x_val = np.array(np.squeeze(ncidin.variables['x0']))
        y_val = np.array(np.squeeze(ncidin.variables['y0']))
        z_val = np.array(np.squeeze(ncidin.variables['z0']))
        dz = z_val[1] - z_val[0]
        lat_val = np.array(np.squeeze(ncidin.variables['lat0']))
        lon_val = np.array(np.squeeze(ncidin.variables['lon0']))
        gm_val = np.array(np.squeeze(ncidin.variables['grid_mapping_0']))
        lon_origin = ncidin.variables['grid_mapping_0'].longitude_of_projection_origin
        lat_origin = ncidin.variables['grid_mapping_0'].latitude_of_projection_origin

        # close radar file
        ncidin.close()
        
        echotop = np.zeros((nt,nmax,ny,nx))+missing_value

        # Find echo tops
        for ii in range(nmax):
            maxht = np.zeros((ny,nx))+missing_value
            for nl in range(nz-1,-1,-1):
                maxht = np.where(np.logical_and((refl[nl,:,:] >= refl_thresh_vals[ii]),
                                                (maxht == missing_value)), nl*dz+dz, maxht)
            echotop[:,ii,:,:] = maxht
            print np.max(maxht)

        # Write echo top values to a file
        ncidout = nc4.Dataset(outFile,'w',format='NETCDF4')
        
        # create dimensions
        time_dim = ncidout.createDimension('time',None) # None implies UNLIMITED
        refl_thresh_dim = ncidout.createDimension('refl_thresh',nmax) 
        y_dim = ncidout.createDimension('y',ny)
        x_dim = ncidout.createDimension('x',nx)

        # create variables
        timeVar = ncidout.createVariable('time',np.float64,('time',),zlib=True )
        xVar = ncidout.createVariable('x',np.float32,('x',),zlib=True )
        yVar = ncidout.createVariable('y',np.float32,('y',),zlib=True )
        latVar = ncidout.createVariable('lat',np.float32,('y','x',),zlib=True )
        lonVar = ncidout.createVariable('lon',np.float32,('y','x',),zlib=True )
        gmVar = ncidout.createVariable('grid_mapping',np.int32)
        refl_thresh = ncidout.createVariable('refl_thresh',np.float32,('refl_thresh',),zlib=True )
        etop = ncidout.createVariable('echo_top',np.float32,('time','refl_thresh','y','x',),
                                      fill_value=missing_value, zlib=True )

        # create variable attributes
        timeVar.standard_name = 'time'
        timeVar.long_name = 'Data time'
        timeVar.units = 'seconds since 1970-01-01T00:00:00Z'
        timeVar.calendar = 'standard'
        timeVar.axis = 'T'
        #timeVar.bounds = 'time_bounds'
        timeVar.comment = datetime
    
        xVar.standard_name = 'projection_x_coordinate'
        xVar.long_name = 'x distance on the projection plane from the origin'
        xVar.units = 'km'
        xVar.axis = 'X'
    
        yVar.standard_name = 'projection_y_coordinate'
        yVar.long_name = 'y distance on the projection plane from the origin'
        yVar.units = 'km'
        yVar.axis = 'Y'
    
        latVar.standard_name = 'latitude'
        latVar.units = 'degrees_north'
    
        lonVar.standard_name = 'longitude'
        lonVar.units = 'degrees_east'
    
        gmVar.grid_mapping_name = 'azimuthal_equidistant'
        gmVar.longitude_of_projection_origin = lon_origin
        gmVar.latitude_of_projection_origin = lat_origin
        gmVar.false_easting = 0
        gmVar.false_northing = 0
    
        # refl_thresh
        refl_thresh.long_name = 'reflectivity_threshold'
        refl_thresh.units = 'dBZ'
        # etop
        etop.long_name = 'reflectivity_maximum_height'
        etop.units = 'km'

        # create global variable attributes
        ncidout.Conventions = Conventions
        ncidout.institution = institution
        ncidout.source = source
        ncidout.title = title
        ncidout.references = references
        ncidout.comment = comment
        currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");
        ncidout.history = 'File created ' + currentTime

        # write vars to file
        timeVar[:] = time_val
        xVar[:] = x_val
        yVar[:] = y_val
        latVar[:] = lat_val
        lonVar[:] = lon_val
        refl_thresh[:] = refl_thresh_vals
        etop[:,:,:,:] = echotop

        ncidout.close()
