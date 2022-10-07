"""
Rain-type code sample wrapper
"""

from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import netCDF4 as nc4
import numpy as np
import os
import raintype as rt
import netcdf_io as net
import logging as log

## ***************** ALGORITHM USER-INPUT PARAMETERS *****************

## reflectivity info
refl_name = 'REFL';
refl_level = 5;
refl_mv = -9999;   #Missing value of reflectivity field.  Only used if not in input file
refl_dx = 1;       #Grid spacing of Cartesian reflectivity data.  Only used if not in input file

## radar info - only use this if data not contained in input file
radar_lat = -0.630447;
radar_lon = 73.10277;

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
truncZconvthres = 43;
dBZformaxconvradius = 46;
weakechothres = 7;
backgrndradius = 5;       #(in km)
maxConvRadius = 10;       #(in km)
minsize = 8;              #(in km^2)
startslope = 50;          #(in km^2)
maxsize = 2000;           #(in km^2)

## Information about where the reflectivity data is located and where outputs should be written.
fileDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf/20111016/';
fileDirOut = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_type_test/20111016/';

## Information about output
title = 'Rain type classification of DYNAMO SPolKa radar data';
institution = 'University of Washington';
source = 'Code used https://github.com/swpowell/raintype_python';
references = 'http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';
comment = 'Based on 2.5km level of interpolated reflectivity data';

## *****************  END USER INPUT PARAMETERS *****************

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for fname in os.listdir(fileDir):
  if fname.endswith('nc'):

    log.info( "file = {}".format(fname) )
  
    #Filename for output
    ncname = str(fileDirOut + 'raintype_' + fname)

    #Open input file
    ncid = nc4.Dataset(str(fileDir + fname),'r')

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
      
    log.info( "outputFormat = {}".format(outputFormat) )
    
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
      missing_value = refl_mv
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
      missing_value = refl_mv

    #Read in reflectivity
    refl = np.array(np.squeeze(ncid.variables[refl_name][:,refl_level-1,:,:]))
    
    #Close input file
    ncid.close()

    #Determine raintype
    (raintype,types) = rt.raintype(refl, refl_missing_val=missing_value, refl_dx=dx,
                                   minZdiff=minZdiff, deepcoszero=deepcoszero,
                                   shallowconvmin=shallowconvmin,truncZconvthres=truncZconvthres,
                                   dBZformaxconvradius=dBZformaxconvradius,
                                   weakechothres=weakechothres, backgrndradius=backgrndradius,
                                   maxConvRadius=maxConvRadius,minsize=minsize,
                                   startslope=startslope, maxsize=maxsize)

    #Output result
    if raintype is not None:
      if outputFormat == 'zeb':
        net.writeZebNetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                           dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                           minsize,startslope,maxsize,title,institution,source,references,
                           comment,bt,toff,lat,lon,alt,dx,dy,dz,raintype,missing_value)
      elif outputFormat == 'cf':
        net.writeCFnetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                          dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                          minsize,startslope,maxsize,title,institution,source,references,
                          comment,tim,x,y,lat,lon,gm,lat_origin,lon_origin,raintype,missing_value)
      else:
        net.writeBasicNetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                             dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                             minsize,startslope,maxsize,title,institution,source,references,
                             comment,dx,radar_lat,radar_lon,raintype.shape[0],raintype.shape[1],
                             raintype,missing_value)
