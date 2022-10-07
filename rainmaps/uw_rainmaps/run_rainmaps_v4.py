""" Raintype and rainmap code wrapper
"""

from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import netCDF4 as nc4
import numpy as np
import os
from uw_raintype import raintype as rt
from uw_raintype import rtfunctions as rtf
from csu_radartools import csu_blended_rain_tropical as cbrt
import rr_functions as rrf

import rr_netcdf_io as io
#import rr_netcdf_io_test as io

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
"""

## ***************** ALGORITHM USER-INPUT PARAMETERS *****************

## reflectivity info
refl_name = 'REFL';
refl_level = 5;
refl_missing_val = -9999;   #Missing value of reflectivity field.  Only used if not in input file
refl_dx = 1;                #Grid spacing of Cartesian reflectivity data.  Only used if not in input file

## other vars required for rainmaps
kdp_name = 'KDP'
kdp_missing_val = -9999
zdr_name = 'ZDR'
zdr_missing_val = -9999

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

## name of raintype variables
rt_min_name = 'rain_type_min'
rt_best_name = 'rain_type_best'
rt_max_name = 'rain_type_max'

## names of rainrate variables
rr_min_name = 'rain_rate_min'
rr_best_name = 'rain_rate_best'
rr_max_name = 'rain_rate_max'
rr_missing = -9999

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

## Create dict of rain rate methods used and RMSE's
## method = 1: R(Kdp, Zdr): RMSE = 1.2 mm/hr
## method = 2: R(Kdp):      RMSE = 2.8 mm/hr
## method = 3: R(Z, Zdr):   RMSE = 1.8 mm/hr
## method = 4: R(Z_all):    RMSE = 5.2 mm/hr
## method = 5: R(Z_c):      RMSE = 5.9 mm/hr
## method = 6: R(Z_s):      RMSE = 1.3 mm/hr
zr_method = {'kdp_zdr':1, 'kdp':2, 'z_zdr':3, 'z_all':4, 'z_conv':5, 'z_stra':6}
rmse = {1:1.2, 2:2.8, 3:1.8, 4:5.2, 5:5.9, 6:1.3}

## rain type input parameters
minZdiff = 20; 
deepcoszero = 40;
shallowconvmin = 28;
truncZconvthres = 42;
dBZformaxconvradius = 45;
weakechothres = 7;
backgrndradius = 5;       #(in km)
maxConvRadius = 10;       #(in km)
minsize = 8;              #(in km^2)
startslope = 50;          #(in km^2)
maxsize = 2000;           #(in km^2)

## Information about where the reflectivity data is located and where outputs should be written.
fileDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf'
fileDirOut = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_rate_sur'

## Output all data (0) or just best data (1)
outputBestOnly = 0

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
#         '20111001','20111002','20111003','20111004','20111005',
#         '20111006','20111007','20111008','20111009','20111010',
#         '20111011','20111012','20111013','20111014','20111015',
#         '20111016','20111017','20111018','20111019','20111020',
#         '20111021','20111022','20111023','20111024','20111025',
#         '20111026','20111027','20111028','20111029','20111030',
#         '20111031',
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
dates = ['20111016','20111123','20120101']

## Information about output
title = 'Rain types and rain maps of DYNAMO SPolKa radar data';
institution = 'University of Washington';
source_uw = 'Code used https://github.com/swpowell/raintype_python';
source_csu = 'Code used https://github.com/CSU-Radarmet/CSU_RadarTools';
references = 'http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';
comment = 'Based on 2.5km level of interpolated radar data';

## *****************  END USER INPUT PARAMETERS *****************

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# derive coefficients for R-Z relationships (R=aaZ^bb)
#aa_conv = np.power(1./a_conv , 1./b_conv)
#bb_conv = 1.0/b_conv
#aa_stra = np.power(1./a_stra , 1./b_stra)
#bb_stra = 1.0/b_stra

# process data for each date
for date in dates:
    print date
    for fname in os.listdir(fileDir+'/'+date):
      if fname.endswith('nc'):
        log.info( "file = {}".format(fname) )
  
        #Filename for output
        outFile = str(fileDirOut+'/'+date+'/rain_' + fname)

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
    
        if not outputBestOnly:

            #------------------
            #CREATE MIN RAINMAP
            #------------------
            # Apply offsets to input data
            refl_min = rtf.ZtoDBZ(rtf.DBZtoZ(refl) - 1)
            kdp_min = kdp - 0.1
            zdr_min = zdr + 0.1

            # Determine raintype
            (rt_min_all,types) = rt.raintype(fname, fileDir+'/'+date, refl_min, missing_value_refl, 
                                             dx, minZdiff, deepcoszero,shallowconvmin,truncZconvthres,
                                             dBZformaxconvradius,weakechothres, backgrndradius,
                                             maxConvRadius,minsize,startslope, maxsize)

            # Create new conv/str/mixed-only rt map (stra=1,3,5, conv=2,4,6) - rt returns with missing values replaced by NaN's
            (rt_min,rt_min_mask) = rrf.assign_conv_stra_mixed(rt_min_all, types['STRATIFORM'], missing_value_refl)
    
            # Run CSU blended rain algorithm and subtract 2*RMSE
            (rr_min,method_min) = cbrt.calc_blended_rain_tropical(refl_min, zdr_min, kdp_min, rt_min, ice_flag=False,
                                                                  fhc=None, predef='True', band='S',
                                                                  thresh_dz=38.0, thresh_zdr=0.25, thresh_kdp=0.3,
                                                                  thresh_frac_ice=0.1, thresh_nexrad=53.0,
                                                                  r_z_a=216., r_z_b=1.39,
                                                                  r_z_a_c=126., r_z_b_c=1.39, r_z_a_s=291., r_z_b_s=1.55,
                                                                  r_kdp_a=59.5202, r_kdp_b=0.7451,
                                                                  r_z_zdr_a=0.0085, r_z_zdr_b=0.9237, r_z_zdr_c=-5.2389/10.,
                                                                  fit_a=None, fit_b=None, method='cr1995',
                                                                  r_kdp_zdr_a=96.5726, r_kdp_zdr_b=0.9315, r_kdp_zdr_c=-2.1140/10.)
    
            # Subtract 2*RSME from non-zero values
            for i in np.arange(1,len(zr_method)+1):
                zero_mask =  np.zeros(rr_min.shape)
                for x in range(0,rr_min.shape[0]):
                    for y in range(0,rr_min.shape[1]):
                        if (method_min[x,y] == i) and (rr_min[x,y] <= (2*rmse[i]) ):
                            zero_mask[x,y] = 1
                rr_min[method_min == i] = rr_min[method_min == i] - (2*rmse[i])
                rr_min[zero_mask==1] = 0

            #for i in np.arange(1,len(zr_method)+1):
            #    zero_mask =  np.zeros(rr_min.shape)
            #    zero_mask[(method_min == i) and (rr_min <= rmse[i])] = 1
            #    rr_min[method_min == i] = rr_min[method_min == i] - rmse[i]
            #    rr_min[zero_mask==1] = 0

            # TESTING
            #f = open('rr_min2.txt','w')
            #for x in range(0,300):
            #    rr_min[x,:].tofile(f,sep='\t',format='%f')
            #    f.write('\n')
            #f.close()

            # Reapply missing value mask to rt_min_all and rr_min
            rt_min_all[rt_min_mask==1] = missing_value_refl
            rr_min[rt_min_mask==1] = missing_value_refl

        #-------------------
        #CREATE BEST RAINMAP
        #-------------------
        # Determine raintype
        (rt_best_all,types) = rt.raintype(fname, fileDir+'/'+date, refl, missing_value_refl, 
                                          dx, minZdiff, deepcoszero,shallowconvmin,truncZconvthres,
                                          dBZformaxconvradius,weakechothres, backgrndradius,
                                          maxConvRadius,minsize,startslope, maxsize)
     
        # Create new conv/str/mixed-only rt map (stra=1,5, conv=2,4,6, mixed=3) - rt returns with missing values replaced by NaN's
        (rt_best,rt_best_mask) = rrf.assign_conv_stra_mixed(rt_best_all, types['UNCERTAIN'], missing_value_refl)
    
        # Run CSU blended rain algorithm
        (rr_best,method_best) = cbrt.calc_blended_rain_tropical(refl, zdr, kdp, rt_best, ice_flag=False,
                                                                fhc=None, predef='True', band='S',
                                                                thresh_dz=38.0, thresh_zdr=0.25, thresh_kdp=0.3,
                                                                thresh_frac_ice=0.1, thresh_nexrad=53.0,
                                                                r_z_a=216., r_z_b=1.39,
                                                                r_z_a_c=126., r_z_b_c=1.39, r_z_a_s=291., r_z_b_s=1.55,
                                                                r_kdp_a=59.5202, r_kdp_b=0.7451,
                                                                r_z_zdr_a=0.0085, r_z_zdr_b=0.9237, r_z_zdr_c=-5.2389/10.,
                                                                fit_a=None, fit_b=None, method='cr1995',
                                                                r_kdp_zdr_a=96.5726, r_kdp_zdr_b=0.9315, r_kdp_zdr_c=-2.1140/10.)
    
        # Reapply missing value mask to rt_best and rr_best
        rt_best_all[rt_mask==1] = missing_value_refl
        rr_best[rt_mask==1] = missing_value_refl
    
        if not outputBestOnly:

            #------------------
            #CREATE MAX RAINMAP
            #------------------
            # Apply offsets to input data
            refl_max = rtf.ZtoDBZ(rtf.DBZtoZ(refl) + 1)
            kdp_max = kdp + 0.1
            zdr_max = zdr - 0.1
    
            # Determine raintype
            (rt_max_all,types) = rt.raintype(fname, fileDir+'/'+date, refl_max, missing_value_refl, 
                                             dx, minZdiff, deepcoszero,shallowconvmin,truncZconvthres,
                                             dBZformaxconvradius,weakechothres, backgrndradius,
                                             maxConvRadius,minsize,startslope, maxsize)

            # Create new conv/str/mixed-only rt map (stra=1,5, conv=2,3,4,6) - rt returns with missing values replaced by NaN's
            (rt_max,rt_mask) = rrf.assign_conv_stra_mixed(rt_max_all, types['CONVECTIVE'], missing_value_refl)

            # Run CSU blended rain algorithm, then apply convective ZR to mixed type
            (rr_max,method_max) = cbrt.calc_blended_rain_tropical(refl_max, zdr_max, kdp_max, rt_max, ice_flag=False,
                                                                  fhc=None, predef='True', band='S',
                                                                  thresh_dz=38.0, thresh_zdr=0.25, thresh_kdp=0.3,
                                                                  thresh_frac_ice=0.1, thresh_nexrad=53.0,
                                                                  r_z_a=216., r_z_b=1.39,
                                                                  r_z_a_c=126., r_z_b_c=1.39, r_z_a_s=291., r_z_b_s=1.55,
                                                                  r_kdp_a=59.5202, r_kdp_b=0.7451,
                                                                  r_z_zdr_a=0.0085, r_z_zdr_b=0.9237, r_z_zdr_c=-5.2389/10.,
                                                                  fit_a=None, fit_b=None, method='cr1995',
                                                                  r_kdp_zdr_a=96.5726, r_kdp_zdr_b=0.9315, r_kdp_zdr_c=-2.1140/10.)

            # Add 2*RSME to non-zero values
            for i in np.arange(1,len(zr_method)+1):
                zero_mask =  np.zeros(rr_max.shape)
                for x in range(0,rr_max.shape[0]):
                    for y in range(0,rr_max.shape[1]):
                        if (method_max[x,y] == i) and (rr_max[x,y] <= 0.0005):
                            zero_mask[x,y] = 1
                rr_max[method_max == i] = rr_max[method_max == i] + rmse[i]
                rr_max[zero_mask==1] = 0

            # Reapply missing value mask to rt_max and rr_max
            rt_max_all[rt_mask==1] = missing_value_refl
            rr_max[rt_mask==1] = missing_value_refl

        #Output results
        if not outputBestOnly:
            if rt_min is not None and rt_best is not None and rt_max is not None:
                if outputFormat == 'zeb':
                    io.writeZebNetcdf(outFile,rt_min_all,rt_best_all,rt_max_all,rr_min,rr_best,rr_max,
                                      method_min,method_best,method_max,
                                      deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                      weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                      a_conv,b_conv,a_stra,b_stra,
                                      title,institution,source_uw,source_csu,references,comment,
                                      bt,toff,lat,lon,alt,dx,dy,dz,types,missing_value_refl)
                elif outputFormat == 'cf':
                    io.writeCFnetcdf(outFile,rt_min_all,rt_best_all,rt_max_all,rr_min,rr_best,rr_max,
                                     method_min,method_best,method_max,
                                     deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                     weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                     a_conv,b_conv,a_stra,b_stra,
                                     title,institution,source_uw,source_csu,references,comment,
                                     tim,x,y,lat,lon,gm,lat_origin,lon_origin,types,missing_value_refl)
                else:
                    io.writeBasicNetcdf(outFile,rt_min_all,rt_best_all,rt_max_all,rr_min,rr_best,rr_max,
                                        method_min,method_best,method_max,
                                        deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                        weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                        title,institution,source_uw,source_csu,references,comment,
                                        dx,radar_lat,radar_lon,raintype.shape[0],rtout.shape[1],types,missing_value_refl)
        else:
            if rt_best is not None:
                if outputFormat == 'zeb':
                    io.writeCFnetcdfBestOnly(outFile,rt_best_all,rr_best,method_best,
                                             deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                             weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                             a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                             a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                             title,institution,source_uw,source_csu,references,comment,
                                             tim,x,y,lat,lon,gm,lat_origin,lon_origin,types,missing_value_refl)
                elif outputFormat == 'cf':
                    io.writeCFnetcdfBestOnly(outFile,rt_best_all,rr_best,method_best,
                                             deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                             weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                             a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                             a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                             title,institution,source_uw,source_csu,references,comment,
                                             tim,x,y,lat,lon,gm,lat_origin,lon_origin,types,missing_value_refl)
                else:
                    io.writeCFnetcdfBestOnly(outFile,rt_best_all,rr_best,method_best,
                                             deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                                             weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                                             a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                                             a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                                             title,institution,source_uw,source_csu,references,comment,
                                             tim,x,y,lat,lon,gm,lat_origin,lon_origin,types,missing_value_refl)
