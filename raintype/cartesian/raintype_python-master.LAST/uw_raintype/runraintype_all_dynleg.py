""" Rain-type code wrapper
"""

from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import netCDF4 as nc4
import numpy as np
import os
from uw_raintype import raintype as rt
from uw_raintype import netcdf_io as net
from uw_raintype import rtfunctions as rtf
import logging as log  

"""
This code assumes that the input reflectivity file contains 'refl' in its name.
The output ncfile will replace the 'refl' with 'raintype' and leave the rest of the 
   filename the same.

The variables listed in the left column immediately below are those in the user-input
  parameters farther below. The variables listed in the right column below are the
  corresponding variable names in Table 1 of Powell et al. (2016).

  Variable name in this code          Variable name in Powell et al.
  --------------------------          ------------------------------
  minZdiff                            a
  deepcoszero                         b
  shallowconvmin                      Z_shallow
  truncZconvthres                     Z_th
  dBZformaxconvradius                 Z_conv
  weakechothres                       Z_weak
  backgrndradius                      R_bg
  maxconvRadius                       R_conv
  minsize                             A_low
  startslope                          A_med
  maxsize                             A_high

  Inputs:
  refl = Reflectivity
  refl_missing_val = missing value in reflectivity data
  refl_dx (km) = horizontal spacing of grid
  minZdiff = factor for comparing echo to background reflectivity; see equation (1) in journal 
     article referenced above
  deepcoszero = see equation (1) in journal article referenced above
  shallowconvmin = minimum dBZ for classification as convective for objects with area less than 
     startslope
  truncZconvthres = reflectivity threshold at or above which echos are classified as convective;  
     The value in Powell et al. (2016) was used for an S-band radar with a beam width of 0.91 degrees. 
     For C-band and/or larger beam width, this value will probably need to be decreased.  Rain
     type classification is most sensitive to this input.
  dBZformaxconvradius = minimum dBZ required for max_conv_radius to apply; should be somewhere close 
     to truncZconvthres
  weakechothres = minimum dBZ for classification as not weak echo; don't change this without a good 
     reason.  7 is about as low as we can go without getting into Bragg scatter territory.
  backgrndradius (km) = radius within which background reflectivity is computed
  maxConvRadius (km) = maximum radius around convective core for possible uncertain classification; 
     Powell et al. (2016) tested 5, and showed that too much convection was included 
     in stratiform region.  Don't lower this number without a good reason.
  minsize (km^2) = minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV
     classification (See dBZcluster for use)
  startslope (km^2) = any contiguous echo object with areal coverage greater than this but less than 
     maxsize gets a new convective threshold that is linearly interpolated between shallowconvmin and 
     truncZconvthres depending on where between startslope and maxsize its area is (See makedBZcluster)
  maxsize (km^2) = any contiguous echo object greater than this size gets a convective threshold of 
     truncZconvthres (See makedBZcluster)
"""

## ***************** ALGORITHM USER-INPUT PARAMETERS *****************

## reflectivity info
refl_name = 'REFL';
refl_level = 5;
refl_missing_val = -9999;   #Missing value of reflectivity field.  Only used if not in input file
refl_dx = 1;       #Grid spacing of Cartesian reflectivity data.  Only used if not in input file

## radar info - only use this if data not contained in input file
radar_lat = -0.630447;
radar_lon = 73.10277;
min_radius = 12.5
max_radius = 147.

## preferred netcdf output format - one of 'basic', 'cf' (CF compliant) or 'zeb' (Zebra compliant)
## NOTE: if the input file does not contain the fields required for the preferred output format
## then the output format will be set to 'basic'; if you are unsure, leave this set to 'cf'
outputFormat = 'cf'

## variables required in input file for cf or zebra compliant output; these names may be spelled
## slightly differently from time to time so we include them as inputs
## ********** DO NOT REMOVE ANY ELEMENTS OR CHANGE THE ORDER OF THE ARRAYS **********
var_cf = ['time','x0','y0','lat0','lon0','grid_mapping_0']
var_zeb = ['base_time','time_offset','lat','lon','alt','x_spacing','y_spacing','z_spacing']

## rain type input parameters
minZdiff = 20; 
deepcoszero = 40;
shallowconvmin = 28;
truncZconvthres = 38;   # previous value was 42
dBZformaxconvradius = truncZconvthres + 5;   # previous value was 45
weakechothres = 7;
backgrndradius = 5;       #(in km)
maxConvRadius = 10;       #(in km)
minsize = 8;              #(in km^2)
startslope = 50;          #(in km^2)
maxsize = 2000;           #(in km^2)

## Information about where the reflectivity data is located and where outputs should be written.
fileDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf/refl'
fileDirOut = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_type_sur'

## Information about output
title = 'Rain type classification of DYNAMO SPolKa radar data';
institution = 'University of Washington';
source = 'Code used https://github.com/swpowell/raintype_python';
references = 'http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';
comment = 'Based on 2.5km level of interpolated reflectivity data';

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

## Information about output
institution = 'University of Washington';
source = 'SPolKa radar data';
title = 'Rain type classifications';
references1 = 'http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';
references2 = 'Code used https://github.com/swpowell/raintype_python';
comment = 'Based on 2.5km level of interpolated reflectivity data';

## *****************  END USER INPUT PARAMETERS *****************

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for date in dates:
    
  print date
  fileDirDate = fileDir+'/'+date
    
  for fname in os.listdir(fileDirDate):
    if fname.endswith('nc'):

      log.info( "file = {}".format(fname) )
  
      #Filename for output
      ncname = str(fileDirOut+'/'+date+'/'+fname.replace('refl','raintype') )
    
      #If output dir does not exist, create it
      dir = os.path.dirname(ncname)
      if not os.path.exists(dir):
        os.makedirs(dir)

      #Open input file
      ncid = nc4.Dataset(str(fileDirDate+'/'+fname),'r')

      #If check to make sure all vars necessary for output format are present in input file
      if outputFormat == 'zeb':
        for x in var_zeb:
          try:
            ncid.variables[x]
          except:
            outputFormat = 'basic'
            break
      elif outputFormat == 'cf':  
        for x in var_cf:
          try:
            ncid.variables[x]
          except:
            outputFormat = 'basic'
            break
      
      #log.info( "outputFormat = {}".format(outputFormat) )
    
      #Make sure refl_name exists.
      try:
        ncid.variables[refl_name]
      except:
        raise SystemExit('Name of reflectivity variable is incorrect. See user input.')
   
      #If variables required by zebra netcdf files are present, read them
      if outputFormat == 'zeb':
        bt = ncid.variables[var_zeb[0]][:]
        toff = ncid.variables[var_zeb[1]][:]
        lat = ncid.variables[var_zeb[2]][:]
        lon = ncid.variables[var_zeb[3]][:]
        alt = ncid.variables[var_zeb[4]][:]
        dx = ncid.variables[var_zeb[5]][:]
        dy = ncid.variables[var_zeb[6]][:]
        dz = ncid.variables[var_zeb[7]][:]
        missing_value = ncid.variables[refl_name].missing_value
      #If variables required for cf-compliancy are present, read them
      elif outputFormat == 'cf':
        dx = refl_dx
        missing_value = refl_missing_val
        tim = ncid.variables[var_cf[0]][:]
        x = ncid.variables[var_cf[1]][:]
        y = ncid.variables[var_cf[2]][:]
        lat = ncid.variables[var_cf[3]][:]
        lon = ncid.variables[var_cf[4]][:]
        gm = ncid.variables[var_cf[5]][:]
        #gmAtts = ncid.variables[var_cf[5]]
        lat_origin = ncid.variables[var_cf[5]].latitude_of_projection_origin
        lon_origin = ncid.variables[var_cf[5]].longitude_of_projection_origin
      else:
        dx = refl_dx
        missing_value = refl_missing_val

      #Read in reflectivity
      refl = np.array(np.squeeze(ncid.variables[refl_name][:,refl_level-1,:,:]))

      #Close input file
      ncid.close()

      #Determine raintype
      (rtout,types) = rt.raintype(fname, fileDirDate, refl, refl_missing_val=missing_value, 
                                  refl_dx=dx, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                  shallowconvmin=shallowconvmin,truncZconvthres=truncZconvthres,
                                  dBZformaxconvradius=dBZformaxconvradius,
                                  weakechothres=weakechothres, backgrndradius=backgrndradius,
                                  maxConvRadius=maxConvRadius,minsize=minsize,
                                  startslope=startslope, maxsize=maxsize)

      # Apply range mask
      (ydim,xdim) = refl.shape
      mask = rtf.radial_distance_mask(min_radius,max_radius,xdim,ydim,dx,dx)
      for jj in range(0,ydim):
        for ii in range(0,xdim):
          if np.logical_and( rtout[jj,ii]==missing_value,mask[jj,ii]==1 ):
            rtout[jj,ii] = 0

      #Output result
      if rtout is not None:
        if outputFormat == 'zeb':
          net.writeZebNetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                             dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                             minsize,startslope,maxsize,title,institution,source,references1,
                             references2,comment,bt,toff,lat,lon,alt,dx,dy,dz,rtout,missing_value)
        elif outputFormat == 'cf':
          net.writeCFnetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                            dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                            minsize,startslope,maxsize,title,institution,source,references1,
                            references2,comment,tim,x,y,lat,lon,gm,lat_origin,lon_origin,rtout,
                            missing_value)
        else:
          net.writeBasicNetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                               dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                               minsize,startslope,maxsize,title,institution,source,references1,
                               references2,comment,dx,radar_lat,radar_lon,raintype.shape[0],
                               rtout.shape[1],rtout,missing_value)
