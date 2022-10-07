""" Rain-type code wrapper
"""

from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import netCDF4 as nc4
import numpy as np
import os
import sys
from rel_raintype import raintype as rt
from rel_raintype import netcdf_io as net
from rel_raintype import rtfunctions as rtf
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

## Read input parameter
if len(sys.argv) != 2:
    print >>sys.stderr, "Useage: ",sys.argv[0]," ['SP'|'MS']"
    quit()
author = sys.argv[1]

## ***************** ALGORITHM USER-INPUT PARAMETERS *****************

## reflectivity info
refl_name = 'DBZ';
refl_level = 4;
refl_missing_val = -9999;   #Missing value of reflectivity field.  Only used if not in input file
refl_dx = 0.75;       #Grid spacing of Cartesian reflectivity data.  Only used if not in input file

## radar info - only use this if data not contained in input file
radar_lat = -31.44133;
radar_lon = -64.19192;
min_radius = 12.5
max_radius = 235.

## preferred netcdf output format - one of 'basic', 'cf' (CF compliant) or 'zeb' (Zebra compliant)
## NOTE: if the input file does not contain the fields required for the preferred output format
## then the output format will be set to 'basic'; if you are unsure, leave this set to 'cf'
outputFormat = 'cf_rel'

## variables required in input file for cf or zebra compliant output; these names may be spelled
## slightly differently from time to time so we include them as inputs
## ********** DO NOT REMOVE ANY ELEMENTS OR CHANGE THE ORDER OF THE ARRAYS **********
var_cf = ['time','x0','y0','lat0','lon0','grid_mapping_0']
var_zeb = ['base_time','time_offset','lat','lon','alt','x_spacing','y_spacing','z_spacing']
var_cf_rel = ['time','x0','y0','z0','grid_mapping_0']

## rain type input parameters
minZdiff = 20
deepcoszero = 40
shallowconvmin = 28
truncZconvthres = 38   # previous value was 42
dBZformaxconvradius = truncZconvthres + 5   # previous value was 45
weakechothres = 7
backgrndradius = 5       #(in km)
if author == 'SP':
  maxConvRadius = 10     #(in km)
elif author == 'MS':
  maxConvRadius = 5      #(in km)
minsize = 8              #(in km^2)
startslope = 50          #(in km^2)
maxsize = 2000           #(in km^2)
print >>sys.stderr, "maxConvRadius = ", maxConvRadius

## Information about where the reflectivity data is located and where outputs should be written.
fileDir = '/home/storm/brodzik/data/radarCart/RMA1/long_range/20181111'
fileDirOut = '/home/storm/brodzik/data/raintype/RMA1/long_range/steiner_comp/20181111'

## Information about output
institution = 'University of Washington';
source = 'Cordoba AR radar data';
title = 'Rain type classifications';
references1 = 'http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';
references2 = 'Code used https://github.com/swpowell/raintype_python';
comment = 'Based on 2.5km level of interpolated reflectivity data';

## *****************  END USER INPUT PARAMETERS *****************

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

log.info( "Algorithm author = {}".format(author) )

for fname in os.listdir(fileDir):
  if fname.endswith('nc'):

    log.info( "file = {}".format(fname) )
  
    #Filename for output
    ncname = str(fileDirOut+'/'+fname.replace('ncf','raintype') )
    
    #If output dir does not exist, create it
    dir = os.path.dirname(ncname)
    if not os.path.exists(dir):
      os.makedirs(dir)

    #Open input file
    ncid = nc4.Dataset(str(fileDir+'/'+fname),'r')

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
    elif outputFormat == 'cf_rel':  
      for x in var_cf_rel:
        try:
          ncid.variables[x]
        except:
          outputFormat = 'basic'
          break
      
    #log.info( "outputFormat = {}".format(outputFormat) )
    
    #Make sure refl_name exists.
    try:
      ncid.variables[refl_name]
    except Exception,e:
      print >>sys.stderr, 'field name', e, 'does not exist'
      continue
   
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
    elif outputFormat == 'cf_rel':
      dx = refl_dx
      missing_value = refl_missing_val
      tim = ncid.variables[var_cf_rel[0]][:]
      x = ncid.variables[var_cf_rel[1]][:]
      y = ncid.variables[var_cf_rel[2]][:]
      z = ncid.variables[var_cf_rel[3]][:]
      #lat = ncid.variables[var_cf_rel[3]][:]
      #lon = ncid.variables[var_cf_rel[4]][:]
      gm = ncid.variables[var_cf_rel[4]][:]
      #gmAtts = ncid.variables[var_cf_rel[3]]
      lat_origin = ncid.variables[var_cf_rel[4]].latitude_of_projection_origin
      lon_origin = ncid.variables[var_cf_rel[4]].longitude_of_projection_origin
    else:
      dx = refl_dx
      missing_value = refl_missing_val

    #Read in reflectivity - ORIGINAL WAY USING 2.5km LEVEL
    #refl_level = 4;  #2.5km
    #refl_level = 7;   #4.0km
    refl = np.array(np.squeeze(ncid.variables[refl_name][:,refl_level-1,:,:]))

    #Read in reflectivity and create 2d array of max in each vertical column
    #refl3d = np.array(np.squeeze(ncid.variables[refl_name][:]))
    #refl = np.max(refl3d,axis=0)
    
    #Read in reflectivity and create 2d array of max in each vertical column
    #refl3d = np.array(np.squeeze(ncid.variables[refl_name][:]))
    #tempArr = np.copy(refl3d)
    #tempArr[tempArr == refl_missing_val] = np.nan
    #counts = np.count_nonzero(~np.isnan(tempArr),axis=0)
    #refl = np.max(refl3d,axis=0)
    #refl[counts<9] = refl_missing_val
    
    #Read in reflectivity below 4.5km and create 2d array of max in each vertical column
    #index range below BB
    #level1 = 0
    #level2 = 7
    #index range above BB
    #level1 = 9
    #level2 = 26
    #ht_below_bb = 4.
    #ht_above_bb = 5.5
    #refl3d = np.array(np.squeeze(ncid.variables[refl_name][:,level1:level2,:,:]))
    #print >>sys.stderr, "refl3d.shape = ", refl3d.shape
    #refl = np.max(refl3d,axis=0)
    
    #Read in reflectivity below 5km and create 2d array of max in each vertical column
    #level = 8
    #refl3d = np.array(np.squeeze(ncid.variables[refl_name][:,0:level,:,:]))
    #tempArr = np.copy(refl3d)
    #tempArr[tempArr == refl_missing_val] = np.nan
    #counts = np.count_nonzero(~np.isnan(tempArr),axis=0)
    #refl = np.max(refl3d,axis=0)
    #refl[counts<3] = refl_missing_val
    
    #Close input file
    ncid.close()

    #Determine raintype
    (rtout,types) = rt.raintype(fname, fileDir, refl, refl_missing_val=missing_value, 
                                refl_dx=dx, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                shallowconvmin=shallowconvmin,truncZconvthres=truncZconvthres,
                                dBZformaxconvradius=dBZformaxconvradius,
                                weakechothres=weakechothres, backgrndradius=backgrndradius,
                                maxConvRadius=maxConvRadius,minsize=minsize,
                                startslope=startslope, maxsize=maxsize)

    #Reclassify extra categories for MS algorith
    #types = {'NO_ECHO':0,'STRATIFORM':1,'CONVECTIVE':2,'UNCERTAIN':3,'ISO_CONV_CORE':4,
    #         'ISO_CONV_FRINGE':5,'WEAK_ECHO':6,'CS_CORE':8,'ISO_CS_CORE':9}
    #log.info( "Before reclassify loop" )
    #log.info( "  rtout values = {}".format(np.unique(rtout)) )
    if author == 'MS':
      rtout[(rtout==3)] = 2
      rtout[(rtout==4)] = 2
      rtout[(rtout==5)] = 1
      rtout[(rtout==6)] = 3
    #log.info( "After reclassify loop" )
    #log.info( "  rtout values = {}".format(np.unique(rtout)) )

    #Output result
    if rtout is not None:
      if outputFormat == 'zeb':
        net.writeZebNetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                           dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                           minsize,startslope,maxsize,title,institution,source,references,
                           comment,bt,toff,lat,lon,alt,dx,dy,dz,rtout,missing_value)
      elif outputFormat == 'cf':
        net.writeCFnetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                          dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                          minsize,startslope,maxsize,title,institution,source,references,
                          comment,tim,x,y,lat,lon,gm,lat_origin,lon_origin,rtout,missing_value)
      elif outputFormat == 'cf_rel':
        net.writeCFRELnetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                             dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                             minsize,startslope,maxsize,title,institution,source,references1,
                             references2,comment,tim,x,y,gm,lat_origin,lon_origin,rtout,
                             missing_value,author)
      else:
        net.writeBasicNetcdf(ncname,types,deepcoszero,shallowconvmin,minZdiff,truncZconvthres,
                             dBZformaxconvradius,weakechothres,backgrndradius,maxConvRadius,
                             minsize,startslope,maxsize,title,institution,source,references,
                             comment,dx,radar_lat,radar_lon,raintype.shape[0],rtout.shape[1],
                             rtout,missing_value)
