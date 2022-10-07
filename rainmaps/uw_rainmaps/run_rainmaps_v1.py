# ***** Rainrate code for DYNAMO legacy products*****
# Author: Stacy Brodzik, University of Washington
# Date: April 26, 2016
# Description: 

from __future__ import absolute_import
from __future__ import division
from csu_radartools import csu_blended_rain_gan as cbrg
import netCDF4 as nc4
import os
import numpy as np
import logging as log
import math
import rr_netcdf_io as net
#from cdf_lib import rr_netcdf_io as net

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
reflDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf/'
rtDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_type_sur/'
rrDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_rate_sur/'

## dates to process
#date = ['20111001','20111002','20111003','20111004','20111005','20111006','20111007',  \
#        '20111008','20111009','20111010','20111011','20111012','20111013','20111014',  \
#        '20111015','20111016','20111017','20111018','20111019','20111020','20111021',  \
#        '20111022','20111023','20111024','20111025','20111026','20111027','20111028',  \
#        '20111029','20111030','20111031','20111101','20111102','20111103','20111104',  \
#        '20111105','20111106','20111107','20111108','20111109','20111110','20111111',  \
#        '20111112','20111113','20111114','20111115','20111116','20111117','20111118',  \
#        '20111119','20111120','20111121','20111122','20111123','20111124','20111125',  \
#        '20111126','20111127','20111128','20111129','20111130','20111201','20111202',  \
#        '20111203','20111204','20111205','20111206','20111207','20111208','20111209',  \
#        '20111210','20111211','20111212','20111213','20111214','20111215','20111216',  \
#        '20111217','20111218','20111219','20111220','20111221','20111222','20111223',  \
#        '20111224','20111225','20111226','20111227','20111228','20111229','20111230',  \
#        '20111231','20120101','20120102','20120103','20120104','20120105','20120106',  \
#        '20120107','20120108','20120109','20120110','20120111','20120112','20120113',  \
#        '20120114','20120115','20120116']
date = ['20111016']

## reflectivity info
refl_name = 'REFL'
refl_level = 5
refl_mv = -9999;   #Missing value of reflectivity field.  Only used if not in input file
refl_dx = 1;       #Grid spacing of Cartesian reflectivity data.  Only used if not in input file

## polarimetric info needed for CSU blended approach
zdr_name = 'ZDR'
kdp_name = 'KDP'
calc_csu_blended = True

## radar info - only use this if data not contained in input file
radar_lat = -0.630447;
radar_lon = 73.10277;

## preferred netcdf output format
#  one of 'basic', 'cf' (CF compliant) or 'zeb' (Zebra compliant)
#  NOTE: if the input file does not contain the fields required for the preferred output format
#  then the output format will be set to 'basic'; if you are unsure, leave this set to 'cf'
outputFormat = 'cf'

## variables required in input file for cf or zebra compliant output; these names may be spelled
## slightly differently from time to time so we include them as inputs
## ********** DO NOT REMOVE ANY ELEMENTS OR CHANGE THE ORDER OF THE ARRAYS **********
var_cf = ['time','x0','y0','lat0','lon0','grid_mapping_0']
var_zeb = ['base_time','time_offset','lat','lon','alt','x_spacing','y_spacing','z_spacing']

# name of raintype variable
rt_name = 'rain_type'

# names of rainrate variables
rr_min_name = 'rain_rate_min'
rr_max_name = 'rain_rate_max'
rr_missing = -9999

#Z-R (Z=aR^b) coefficients for MISMO single equation estimation (Z in linear units mm^6/m^3)
a_mismo = 178.0
b_mismo = 1.44

# Z-R coefficients and exponents for dual equation estimation
a_conv = 134.0
b_conv = 1.44
a_stra = 300.28
b_stra = 1.55

## Information about output
title = 'Rain rate maps based on DYNAMO SPolKa radar data';
institution = 'University of Washington';
source_uw = 'Code used https://github.com/swpowell/rainrate_python';
source_csu = 'Code used https://github.com/CSU-Radarmet/CSU_RadarTools';
references = 'http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';
comment = 'Based on 2.5km level of interpolated reflectivity data. Min rate includes UNCERTAIN raintype as stratiform.  Max rate includes UNCERTAIN raintype as convective';

# set up logging
log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)
# ---------------------------------- DONE WITH INPUTS ----------------------------------

# derive coefficients for R-Z relationships (R=aaZ^bb)
aa_mismo = np.power(1./a_mismo , 1./b_mismo)
bb_mismo = 1.0/b_mismo
aa_conv = np.power(1./a_conv , 1./b_conv)
bb_conv = 1.0/b_conv
aa_stra = np.power(1./a_stra , 1./b_stra)
bb_stra = 1.0/b_stra

# process data for each date
for i,idate in enumerate(date):
    print idate
    for reflFile in os.listdir(reflDir+idate):
        if reflFile.endswith('nc'):
            log.info('reflFile = {}'.format(reflFile) )

            # open input reflectivity file
            ncid_refl = nc4.Dataset(str(reflDir + idate + '/' + reflFile), 'r')

            # assign output rainrate file name
            rrFile = str(rrDir + idate + '/rainrate_' + reflFile)

            #If output dir does not exist, create it
            dir = os.path.dirname(rrFile)
            if not os.path.exists(dir):
                os.makedirs(dir)

            # check to make sure all vars necessary for output format are present in input file
            if outputFormat == 'zeb':
                for x in var_zeb:
                    try:
                        ncid_refl.variables[x]
                    except:
                        outputFormat = 'basic'
                        break
            elif outputFormat == 'cf':  
                for x in var_cf:
                    try:
                        ncid_refl.variables[x]
                    except:
                        outputFormat = 'basic'
                        break
      
            log.info( "   outputFormat = {}".format(outputFormat) )
    
            # check for existence of reflectivity
            try:
                ncid_refl.variables[refl_name]
                log.info( "   {} present".format(refl_name) )
            except:
                raise SystemExit('Name of reflectivity variable is incorrect or missing.')

            # check for existence of ZDR and KDP
            for x in [zdr_name,kdp_name]:
                try:
                    ncid_refl.variables[x]
                    log.info( "   {} present".format(x) )
                except:
                    calc_csu_blended = False
                    log.info( "   ZDR and/or KDP missing.  Will not calculate CSU blended rain." )

            log.info( "   calc_csu_blended = {}".format(calc_csu_blended) )    

            # if variables required by zebra netcdf files are present, read them
            if outputFormat == 'zeb':
                bt = ncid_refl.variables[var_zeb[0]][:]
                toff = ncid_refl.variables[var_zeb[1]][:]
                lat = ncid_refl.variables[var_zeb[2]][:]
                lon = ncid_refl.variables[var_zeb[3]][:]
                alt = ncid_refl.variables[var_zeb[4]][:]
                dx = ncid_refl.variables[var_zeb[5]][:]
                dy = ncid_refl.variables[var_zeb[6]][:]
                dz = ncid_refl.variables[var_zeb[7]][:]
                missing_value = ncid_refl.variables[refl_name].missing_value
            #If variables required for cf-compliancy are present, read them
            elif outputFormat == 'cf':
                dx = refl_dx
                missing_value = refl_mv
                tim = ncid_refl.variables[var_cf[0]][:]
                x = ncid_refl.variables[var_cf[1]][:]
                y = ncid_refl.variables[var_cf[2]][:]
                lat = ncid_refl.variables[var_cf[3]][:]
                lon = ncid_refl.variables[var_cf[4]][:]
                gm = ncid_refl.variables[var_cf[5]][:]
                #gmAtts = ncid_refl.variables[var_cf[5]]
                lat_origin = ncid_refl.variables[var_cf[5]].latitude_of_projection_origin
                lon_origin = ncid_refl.variables[var_cf[5]].longitude_of_projection_origin
            else:
                dx = refl_dx
                missing_value = refl_mv

            # read appropriate level of reflectivity data (and polarimetric fields if appropriate)
            refl = np.array(np.squeeze(ncid_refl.variables[refl_name][:,refl_level-1,:,:]))
            if calc_csu_blended:
                zdr = np.array(np.squeeze(ncid_refl.variables[zdr_name][:,refl_level-1,:,:]))
                kdp = np.array(np.squeeze(ncid_refl.variables[kdp_name][:,refl_level-1,:,:]))

            # close reflectivity file
            ncid_refl.close()
            
            # create missing_value mask and turn refl,zdr,kdp missing values into NaN's
            mask = np.zeros(refl.shape)
            mask[(refl == missing_value)] = 1
            refl[(mask == 1)] = np.nan
            if calc_csu_blended:
                zdr[(mask == 1)] = np.nan
                kdp[(mask == 1)] = np.nan

            #Convert dBZ to Z
            Z = 10**(0.1*refl)
            Z_1d = Z.reshape(Z.shape[0]*Z.shape[1])

            # open raintype file -- assume if has same name as refl file with 'raintype_' prepended
            ncid_rt = nc4.Dataset(str(rtDir + idate + '/raintype_' + reflFile), 'r')
            log.info('   rtFile = {}'.format('raintype_'+reflFile) )
            
            # check for existence of raintype
            try:
                ncid_rt.variables[rt_name]
                log.info( "   {} present".format(rt_name) )
            except:
                raise SystemExit('Name of raintype variable is incorrect or missing.')

            # read raintype data
            rt = np.squeeze(ncid_rt.variables[rt_name])

            # close raintype file
            ncid_rt.close()

            # make sure size of refl and size of rt are the same before proceeding
            if refl.shape == rt.shape:

                # apply Z-R relations
                # convection: Z=a_convR^b_conv or R=aa_convZ^bb_conv
                # stratiform: Z=a_straR^b_stra or R=aa_convZ^bb_conv
                # both: Z=a_mismoR^b_mismo or R=aa_mismoZ^bb_mismo

                # min rainmap - apply conv eqn to CONVECTIVE (2), ISO_CONV_CORE (4), ISO_CONV_FRINGE (5)
                #               apply stra eqn to STRATIFORM (1), UNCERTAIN (3)
            
                # max rainmap - apply conv eqn to CONVECTIVE (2), ,UNCERTAIN (3), ISO_CONV_CORE (4),
                #                                 ISO_CONV_FRINGE (5)
                #               apply stra eqn to STRATIFORM (1)

                # raintype breakdown making UNCERTAIN = stratiform
                rt_conv_sU = [2, 4, 5]
                rt_stra_sU = [1, 3]
                # create associated conv/stra masks
                mask_conv_sU = np.in1d(rt,rt_conv_sU)
                mask_stra_sU = np.in1d(rt,rt_stra_sU)
                
                # raintype breakdown making UNCERTAIN = convective
                rt_conv_cU = [2, 3, 4, 5]
                rt_stra_cU = [1]
                # create associated conv/stra masks
                mask_conv_cU = np.in1d(rt,rt_conv_cU)
                mask_stra_cU = np.in1d(rt,rt_stra_cU)

                # calculate 2zr min rainmap
                rr_2zr_min_1d = np.zeros(rt.shape[0]*rt.shape[1])
                rr_2zr_min_1d[mask_conv_sU] = aa_conv * np.power(Z_1d[mask_conv_sU] , bb_conv )
                rr_2zr_min_1d[mask_stra_sU] = aa_stra * np.power(Z_1d[mask_stra_sU] , bb_stra )
                rr_2zr_min = rr_2zr_min_1d.reshape(refl.shape)

                # calculate 2zr max rainmap
                rr_2zr_max_1d = np.zeros(rt.shape[0]*rt.shape[1])
                rr_2zr_max_1d[mask_conv_cU] = aa_conv * np.power(Z_1d[mask_conv_cU] , bb_conv )
                rr_2zr_max_1d[mask_stra_cU] = aa_conv * np.power(Z_1d[mask_stra_cU] , bb_stra )
                rr_2zr_max = rr_2zr_max_1d.reshape(refl.shape)

                # calculate 1zr (MISMO) min rainmap
                rr_1zr_min_1d = np.zeros(rt.shape[0]*rt.shape[1])
                rr_1zr_min_1d[mask_conv_sU] = aa_mismo * np.power(Z_1d[mask_conv_sU] , bb_mismo )
                rr_1zr_min = rr_1zr_min_1d.reshape(refl.shape)

                # calculate 1zr (MISMO) max rainmap
                rr_1zr_max_1d = np.zeros(rt.shape[0]*rt.shape[1])
                rr_1zr_max_1d[mask_conv_cU] = aa_mismo * np.power(Z_1d[mask_conv_cU] , bb_mismo )
                rr_1zr_max = rr_1zr_max_1d.reshape(refl.shape)

                # calculate CSU-blended rainmap - original coefficients
                (rr_csu,method_csu) = cbrg.calc_blended_rain_tropical(refl, zdr, kdp, ice_flag=False,
                                                                      thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3,
                                                                      thresh_frac_ice=0.1, thresh_nexrad=53.0,
                                                                      a=300.28, b=1.55,
                                                                      nex_a=300.28, nex_b=1.55,
                                                                      kdp_zdr_a=94.1767, kdp_zdr_b=0.8765, kdp_zdr_c=-2.0245,
                                                                      dz_zdr_a=0.0124, dz_zdr_b=0.8639, dz_zdr_c=-3.9623,
                                                                      kdp_a=58.9141, kdp_b=0.7522,
                                                                      fit_a=None, fit_b=None, method='cr1995')

                # apply mask
                rr_csu[(mask == 1)] = missing_value
                method_csu[(mask == 1)] = missing_value
                rr_1zr_min[(mask == 1)] = missing_value
                rr_1zr_max[(mask == 1)] = missing_value
                rr_2zr_min[(mask == 1)] = missing_value
                rr_2zr_max[(mask == 1)] = missing_value
                
                # find min and max rr at each grid point
                rr_min = np.minimum(rr_2zr_min,rr_2zr_max)
                rr_min = np.minimum(rr_min,rr_1zr_min)
                rr_min = np.minimum(rr_min,rr_1zr_max)
                rr_min[(mask == 1)] = missing_value

                rr_max = np.maximum(rr_2zr_min,rr_2zr_max)
                rr_max = np.maximum(rr_max,rr_1zr_min)
                rr_max = np.maximum(rr_max,rr_1zr_max)
                rr_max[(mask == 1)] = missing_value

                log.info('   Done calculating rainmaps')
                
                # write rainmaps to netcdf file
                if outputFormat == 'zeb':
                    net.writeZeb_RRNetcdf(rrFile,a_conv,b_conv,a_stra,b_stra,title,
                                          institution,source,references,comment,btVal,
                                          toVal,latVal,lonVal,altVal,xspVal,yspVal,
                                          zspVal,rr_min,rr_max,rr_2zr_min,rr_2zr_max,
                                          rr_1zr_min,rr_1zr_max,rr_csu,method_csu,
                                          missing_value)
                    log.info('   Done writing zeb file')
                elif outputFormat == 'cf':
                    net.writeCF_RRnetcdf(rrFile,a_conv,b_conv,a_stra,b_stra,title,
                                         institution,source,references,comment,
                                         tim,x,y,lat,lon,gm,lat_origin,lon_origin,
                                         rr_min,rr_max,rr_2zr_min,rr_2zr_max,
                                         rr_1zr_min,rr_1zr_max,rr_csu,method_csu,
                                         missing_value)
                    log.info('   Done writing CF file')
                else:
                    net_writeBasic_RRNetcdf(rrFile,a_conv,b_conv,a_stra,b_stra,title,
                                            institution,source,references,comment,dx,
                                            radar_lat,radar_lon,xdim,ydim,rr_min,
                                            rr_max,rr_2zr_min,rr_2zr_max,rr_1zr_min,
                                            rr_1zr_max,rr_csu,method_csu,missing_value)
                    log.info('   Done writing basic file')
                    

            

        
        

    
