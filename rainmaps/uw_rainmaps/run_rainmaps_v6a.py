""" Raintype and rainmap code wrapper
"""

from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import netCDF4 as nc4
import numpy as np
import os
import logging as log
from uw_raintype import raintype as rt
from uw_raintype import rtfunctions as rtf
#from uw_rainmaps import rainmaps as rm
import rainmaps as rm
import rr_functions as rrf
import rr_netcdf_io as io

from copy import copy, deepcopy
#y = deepcopy(x)    to actually copy an array and not just make a reference

import logging as log  

"""
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

Raintypes defined by Powell et al
0 = NO_SFC_ECHO
1 = STRATIFORM
2 = CONVECTIVE
3 = UNCERTAIN
4 = ISO_CONV_CORE
5 = ISO_CONV_FRINGE
6 = WEAK_ECHO
"""

## ***************** ALGORITHM USER-INPUT PARAMETERS *****************

## reflectivity info
refl_name = 'DBZ_S';
refl_level = 5;
refl_missing_val = -9999;   #Missing value of reflectivity field.  Only used if not in input file
refl_dx = 1;                #Grid spacing of Cartesian reflectivity data.  Only used if not in input file

## other vars required for rainmaps
kdp_name = 'KDP_S'
kdp_missing_val = -9999
zdr_name = 'ZDR_S'
zdr_missing_val = -9999

## radar info - only use this if data not contained in input file
radar_lat = -0.630447;
radar_lon = 73.10277;
min_radius = 12.5
max_radius = 147.

## preferred netcdf output format - one of 'basic', 'cf' (CF compliant) or 'zeb' (Zebra compliant)
## NOTE: if the input file does not contain the fields required for the preferred output format
## then the output format will be set to 'basic'; if you are unsure, leave this set to 'cf'
outputFormat = 'zeb'

## variables required in input file for cf or zebra compliant output; these names may be spelled
## slightly differently from time to time so we include them as inputs
## ********** DO NOT REMOVE ANY ELEMENTS OR CHANGE THE ORDER OF THE ARRAYS **********
var_cf = ['time','x0','y0','lat0','lon0','grid_mapping_0']
var_zeb = ['base_time','time_offset','lat','lon','alt','x_spacing','y_spacing','z_spacing']

## names of missing variables
missing_value_rt = -99
missing_value_rr = -99

## Z-R coefficients and exponents for dual equation estimation
a_conv = 126.0
b_conv = 1.46
a_stra = 291.0
b_stra = 1.55
a_all = 216.0
b_all = 1.39
a_kdp_zdr = 96.57
b_KDP_zdr = 0.93
b_kdp_ZDR = -2.11
a_zh_zdr = 0.0085
b_ZH_zdr = 0.92
b_zh_ZDR = -5.24
a_kdp = 56.04
b_kdp = 0.80

## Create dicts of rain rate methods used and sigmaOverR and RMSE values
zr_method = {'kdp_zdr':1, 'kdp':2, 'z_zdr':3, 'z_all':4, 'z_conv':5, 'z_stra':6}
sigOverR = {3:0.307, 4:0.144, 5:0.137, 6:0.129}

uncert_vars = {'a':1, 'b':2, 'sigma_z':3, 'sigma_zdr':4, 'sigma_kdp':5}
uncert_vals = {1:{1:0.932, 2:-2.114, 3:0, 4:0.2, 5:0.8},
               2:{1:0.825, 2:0,      3:0, 4:0,   5:0.8} }

rr_bin = {'lt_20':1, '20-60':2, 'gt_60':3}
rr_bin_stra = {'lt_10':1, '10-20':2, 'gt_20':3}
rmse_a = {1:{1:0.73, 2:0.77, 3:0.94},
          2:{1:0.88, 2:0.63, 3:0.75},
          3:{1:0.32, 2:0.12, 3:0.09},
          4:{1:1.19, 2:0.72, 3:0.95},
          5:{1:0.49, 2:0.21, 3:0.30},
          6:{1:0.78, 2:0.82, 3:0.76} }
rmse_b = {1:{1:0.38, 2:0.37, 3:0.32},
          2:{1:0.57, 2:0.70, 3:0.67},
          3:{1:0.66, 2:0.97, 3:1.06},
          4:{1:0.65, 2:0.83, 3:0.78},
          5:{1:0.80, 2:1.08, 3:1.00},
          6:{1:0.62, 2:0.68, 3:0.78} }

## rain type input parameters
minZdiff = 20; 
deepcoszero = 40;
shallowconvmin = 28;
#truncZconvthres = 42;
truncZconvthres = 38;
#dBZformaxconvradius = 45;
dBZformaxconvradius = truncZconvthres + 5;
weakechothres = 7;
backgrndradius = 5;       #(in km)
maxConvRadius = 10;       #(in km)
minsize = 8;              #(in km^2)
startslope = 50;          #(in km^2)
maxsize = 2000;           #(in km^2)

## Information about where the reflectivity data is located and where outputs should be written.
#fileDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf'
fileDir = '/home/disk/mjo/dynamo/data.server/zebra/QCed/spolka/sur'
#fileDirOut = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_rate_sur/38'
fileDirOut = '/home/disk/mjo/dynamo/data.server/zebra/QCed/spolka/rain_type/sur_0.5km'

## Output all data (0) or just best data (1)
outputBestOnly = 1

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
dates = ['20111026']

## Information about output
institution = 'University of Washington';
source = 'SPolKa radar data';
source_uw = 'Code used https://github.com/swpowell/raintype_python';
source_csu='Code used https://github.com/CSU-Radarmet/CSU_RadarTools';
#title = 'Rainrates: minimum, expected, maximum';
title = 'Rain types and rain maps of DYNAMO SPolKa radar data';
#references = 'Code used https://github.com/CSU-Radarmet/CSU_RadarTools/csu_blended_rain_tropical.py';
references='http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';
comment = 'Based on 2.5km level of interpolated radar data';

## *****************  END USER INPUT PARAMETERS *****************

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# process data for each date
for date in dates:
    print date
    for fname in os.listdir(fileDir+'/'+date):
      if fname.endswith('nc'):
        log.info( "file = {}".format(fname) )
  
        #Filename for output
        outFile = str(fileDirOut+'/'+date+'/rain_' + fname)
        #outFile = str(fileDirOut+'/'+date+'/'+fname.replace('3dvol','rainmap') )

        #Open input file
        ncid = nc4.Dataset(str(fileDir+'/'+date+'/'+fname),'r')
        
        #If output dir does not exist, create it
        dir = os.path.dirname(outFile)
        if not os.path.exists(dir):
          os.makedirs(dir)

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
    
        #Make sure refl, kdp and zdr exist in input file
        for x in [refl_name,zdr_name,kdp_name]:
          try:
            ncid.variables[x]
            log.info( "   {} present".format(x) )
          except:
            log.info( "   refl, ZDR and/or KDP missing.  Exiting." )
            exit()
      
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
          missing_value_refl = ncid.variables[refl_name].missing_value
          missing_value_kdp = ncid.variables[kdp_name].missing_value
          missing_value_zdr = ncid.variables[zdr_name].missing_value
        #If variables required for cf-compliancy are present, read them
        elif outputFormat == 'cf':
          dx = refl_dx
          missing_value_refl = refl_missing_val
          missing_value_kdp = kdp_missing_val
          missing_value_zdr = zdr_missing_val
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
          missing_value_refl = refl_missing_val
          missing_value_kdp = kdp_missing_val
          missing_value_zdr = zdr_missing_val

        #Read in appropriate level of refl, kdp and zdr
        refl = np.array(np.squeeze(ncid.variables[refl_name][:,refl_level-1,:,:]))
        kdp = np.array(np.squeeze(ncid.variables[kdp_name][:,refl_level-1,:,:]))
        zdr = np.array(np.squeeze(ncid.variables[zdr_name][:,refl_level-1,:,:]))

        #Close input file
        ncid.close()

        # create missing_value mask and turn refl,zdr,kdp missing values into NaN's
        refl[(refl==missing_value_refl)] = np.nan
        kdp[(kdp==missing_value_kdp)] = np.nan
        zdr[(zdr==missing_value_zdr)] = np.nan
    
        # Determine raintype
        (rt_powell,types) = rt.raintype(fname, fileDir+'/'+date, refl, missing_value_refl, 
                                        dx, minZdiff, deepcoszero,shallowconvmin,truncZconvthres,
                                        dBZformaxconvradius,weakechothres, backgrndradius,
                                        maxConvRadius,minsize,startslope, maxsize)
        
        # Apply range mask to raintype so that all in-range missing values are set to 0 (no echo)
        (ydim,xdim) = refl.shape
        range_mask = rtf.radial_distance_mask(min_radius,max_radius,xdim,ydim,dx,dx)
        rt_zeros_mask = np.zeros(rt_powell.shape,dtype=np.int)
        for jj in range(0,ydim):
            for ii in range(0,xdim):
                if np.logical_and( rt_powell[jj,ii]==missing_value_refl,range_mask[jj,ii]==1 ):
                    rt_powell[jj,ii] = 0
                    rt_zeros_mask[jj,ii] = 1

        # Determine rainmaps
        (rr_min, rr_best, rr_max, method_min, method_best, method_max) = rm.rainmaps(refl, zdr,
                                        kdp, rt_powell, types, zr_method, sigOverR, uncert_vars, 
                                        uncert_vals, rr_bin, rr_bin_stra, rmse_a, rmse_b,
                                        a_conv, b_conv, a_stra, b_stra, a_all, b_all,
                                        a_kdp_zdr, b_KDP_zdr, b_kdp_ZDR, a_zh_zdr, b_ZH_zdr, b_zh_ZDR,
                                        a_kdp, b_kdp, missing_value_refl, missing_value_rr,
                                        outputBestOnly)
        
        # Ensure that rrate = 0.0 where rtype = 0
        rr_min[rt_zeros_mask==1] = 0.0
        rr_best[rt_zeros_mask==1] = 0.0
        rr_max[rt_zeros_mask==1] = 0.0
        
        #Output results
        #io.writeCFnetcdfV6_rr_only(outFile,rr_min,rr_best,rr_max,
        #                           a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
        #                           a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
        #                           title,institution,source,references,comment,
        #                           tim,x,y,lat,lon,gm,lat_origin,lon_origin,types,missing_value_rr)

        #"""
        if not outputBestOnly:
            if rr_min is not None and rr_best is not None and rr_max is not None:
                if outputFormat == 'zeb':
                    io.writeZebNetcdfV6(outFile,rt_powell,rr_min,rr_best,rr_max,
                                        method_min,method_best,method_max,
                                        deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                        weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                        a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                        a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                        title,institution,source_uw,source_csu,references,comment,
                                        bt,toff,lat,lon,alt,dx,dy,dz,types,missing_value_refl)
                elif outputFormat == 'cf':
                    io.writeCFnetcdfV6(outFile,rt_powell,rr_min,rr_best,rr_max,
                                       method_min,method_best,method_max,zr_method,
                                       deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                       weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                       a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                       a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                       title,institution,source_uw,source_csu,references,comment,
                                       tim,x,y,lat,lon,gm,lat_origin,lon_origin,types,missing_value_refl)
                else:
                    io.writeBasicNetcdfV6(outFile,rt_powell,rr_min,rr_best,rr_max,
                                          method_min,method_best,method_max,
                                          deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                          weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                          a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                          a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                          title,institution,source_uw,source_csu,references,comment,
                                          dx,radar_lat,radar_lon,raintype.shape[0],rtout.shape[1],types,missing_value_refl)
        else:
            if rr_best is not None:
                if outputFormat == 'zeb':
                    io.writeZebNetcdfBestOnlyV6(outFile,rt_powell,rr_best,method_best,
                                                deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                                weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                                a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                                a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                                title,institution,source_uw,source_csu,references,comment,
                                                bt,toff,lat,lon,alt,dx,dy,dz,types,missing_value_refl)
                elif outputFormat == 'cf':
                    io.writeCFnetcdfBestOnlyV6(outFile,rt_powell,rr_best,method_best,
                                               deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                               weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                               a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                               a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                               title,institution,source_uw,source_csu,references,comment,
                                               tim,x,y,lat,lon,gm,lat_origin,lon_origin,types,missing_value_refl)
                else:
                    io.writeBasicNetcdfBestOnlyV6(outFile,rt_powell,rr_best,method_best,
                                                  deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                                  weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                                  a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                                  a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                                  title,institution,source_uw,source_csu,references,comment,
                                                  dx,radar_lat,radar_lon,raintype.shape[0],rtout.shape[1],types,missing_value_refl)
        #"""
