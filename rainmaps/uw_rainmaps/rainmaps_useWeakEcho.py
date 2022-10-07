"""
#Rain-rate code of Brodzik et al.
#Author: Stacy Brodzik
#Date: 1 Feb 2017

#Amendments by

"""

from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import numpy as np
import os
import logging as log
import rr_functions as rrf
from csu_radartools import csu_blended_rain_tropical as cbrt

def rainmaps(refl, zdr, kdp, rt_powell, types, zr_method, sigOverR, uncert_vars, uncert_vals, 
             rr_bin, rr_bin_stra, rmse_a, rmse_b, a_conv, b_conv, a_stra, b_stra, a_all, b_all,
             a_kdp_zdr, b_KDP_zdr, b_kdp_ZDR, a_zh_zdr, b_ZH_zdr, b_zh_ZDR, a_kdp, b_kdp,
             missing_value, missing_value_rr, outputBestOnly):

   Ice_flag=False
   Fhc=None
   Predef='True'
   Band='S'
   Thresh_dz=38.0
   Thresh_zdr=0.25
   Thresh_kdp=0.3
   Thresh_frac_ice=0.1
   Thresh_nexrad=53.0
   Fit_a=None
   Fit_b=None
   Method='cr1995',
   
   if not outputBestOnly:

      #------------------
      #CREATE MIN RAINMAP
      #------------------
      # Create new conv/str/mixed-only rt map (stra=1,3,5, conv=2,4, mixed=6)
      # rt returns with missing values = missing_value (-9999)
      (rt_min,rt_min_mask) = rrf.assign_conv_stra_mixed(rt_powell, types['STRATIFORM'], missing_value)
    
      # Run CSU blended rain algorithm
      (rr_min,method_min) = cbrt.calc_blended_rain_tropical(refl, zdr, kdp, rt_min, ice_flag=Ice_flag,
                                                            fhc=Fhc, predef=Predef, band=Band,
                                                            thresh_dz=Thresh_dz, thresh_zdr=Thresh_zdr,
                                                            thresh_kdp=Thresh_kdp,thresh_frac_ice=Thresh_frac_ice,
                                                            thresh_nexrad=Thresh_nexrad,
                                                            r_z_a=a_all, r_z_b=b_all,
                                                            r_z_a_c=a_conv, r_z_b_c=b_conv,
                                                            r_z_a_s=a_stra, r_z_b_s=b_stra,
                                                            r_kdp_a=a_kdp, r_kdp_b=b_kdp,
                                                            r_z_zdr_a=a_zh_zdr, r_z_zdr_b=b_ZH_zdr, r_z_zdr_c=b_zh_ZDR/10.,
                                                            fit_a=Fit_a, fit_b=Fit_b, method=Method,
                                                            r_kdp_zdr_a=a_kdp_zdr, r_kdp_zdr_b=b_KDP_zdr, r_kdp_zdr_c=b_kdp_ZDR/10.)

      # Subtract 2*RSME from non-zero values
      rmse = np.zeros(rr_min.shape)
      for imethod in np.arange(1,len(zr_method)+1):
         if imethod == 6:
            mask = (method_min == imethod) & (rr_min < 10.)
            bin = rr_bin_stra['lt_10']
            rmse[mask] = rmse_a[imethod][bin] * (rr_min[mask] ** rmse_b[imethod][bin])
                    
            mask = (method_min == imethod) & (rr_min > 20.)
            bin = rr_bin_stra['gt_20']
            rmse[mask] = rmse_a[imethod][bin] * (rr_min[mask] ** rmse_b[imethod][bin])
                    
            mask = (method_min == imethod) & (rr_min >= 10.) & (rr_min <= 20.) 
            bin = rr_bin_stra['10-20']
            rmse[mask] = rmse_a[imethod][bin] * (rr_min[mask] ** rmse_b[imethod][bin])
         else:
            mask = (method_min == imethod) & (rr_min < 20.)
            bin = rr_bin['lt_20']
            rmse[mask] = rmse_a[imethod][bin] * (rr_min[mask] ** rmse_b[imethod][bin])
                    
            mask = (method_min == imethod) & (rr_min > 60.)
            bin = rr_bin['gt_60']
            rmse[mask] = rmse_a[imethod][bin] * (rr_min[mask] ** rmse_b[imethod][bin])
                    
            mask = (method_min == imethod) & (rr_min >= 20.) & (rr_min <= 60.) 
            bin = rr_bin['20-60']
            rmse[mask] = rmse_a[imethod][bin] * (rr_min[mask] ** rmse_b[imethod][bin])
      rr_min = rr_min - (2*rmse)
      rr_min[rr_min < 0.0] = 0.0

      # Subtract standard deviation of measurement error (sigma(epsilon_m)) from non-zero values
      sigmaOverR = np.zeros(rr_min.shape)
      for imethod in np.arange(1,len(zr_method)+1):
         mask = (method_min == imethod)
         if imethod==3 or imethod==4 or imethod==5 or imethod==6:
            sigmaOverR[mask] = sigOverR[imethod]
         else:
            a = uncert_vals[imethod][uncert_vars['a']]
            b = uncert_vals[imethod][uncert_vars['b']]
            sigma_z = uncert_vals[imethod][uncert_vars['sigma_z']]
            sigma_zdr = uncert_vals[imethod][uncert_vars['sigma_zdr']]
            sigma_kdp = uncert_vals[imethod][uncert_vars['sigma_kdp']]
            if imethod==1:
               sigmaOverR[mask] = np.sqrt( (a**2) * (sigma_kdp**2)/(kdp[mask]**2) + 
                                           (b**2) * (sigma_zdr**2)/(zdr[mask]**2) )
            elif imethod==2:
               sigmaOverR[mask] = (a*sigma_kdp) * (a_kdp/rr_min[mask])**(1/a)
      rr_min = rr_min - (sigmaOverR * rr_min)
      rr_min[rr_min <= 0.0] = 0.0
                    
   #-------------------
   #CREATE BEST RAINMAP
   #-------------------
   # Create new conv/str/mixed-only rt map (stra=1,5, conv=2,4, mixed=3,6) - rt returns with missing values = missing_value (-9999)
   (rt_best,rt_best_mask) = rrf.assign_conv_stra_mixed(rt_powell, types['UNCERTAIN'], missing_value)
    
   # Run CSU blended rain algorithm
   (rr_best,method_best) = cbrt.calc_blended_rain_tropical(refl, zdr, kdp, rt_best, ice_flag=Ice_flag,
                                                            fhc=Fhc, predef=Predef, band=Band,
                                                            thresh_dz=Thresh_dz, thresh_zdr=Thresh_zdr,
                                                            thresh_kdp=Thresh_kdp,thresh_frac_ice=Thresh_frac_ice,
                                                            thresh_nexrad=Thresh_nexrad,
                                                            r_z_a=a_all, r_z_b=b_all,
                                                            r_z_a_c=a_conv, r_z_b_c=b_conv,
                                                            r_z_a_s=a_stra, r_z_b_s=b_stra,
                                                            r_kdp_a=a_kdp, r_kdp_b=b_kdp,
                                                            r_z_zdr_a=a_zh_zdr, r_z_zdr_b=b_ZH_zdr, r_z_zdr_c=b_zh_ZDR/10.,
                                                            fit_a=Fit_a, fit_b=Fit_b, method=Method,
                                                            r_kdp_zdr_a=a_kdp_zdr, r_kdp_zdr_b=b_KDP_zdr, r_kdp_zdr_c=b_kdp_ZDR/10.)
   rr_best[rr_best <= 0.0] = 0.0

    
   if not outputBestOnly:

      #------------------
      #CREATE MAX RAINMAP
      #------------------
      # Create new conv/str/mixed-only rt map (stra=1,5, conv=2,3,4, mixed=6) - rt returns with missing values = missing_value (-9999)
      (rt_max,rt_max_mask) = rrf.assign_conv_stra_mixed(rt_powell, types['CONVECTIVE'], missing_value)

      # Run CSU blended rain algorithm, then apply convective ZR to mixed type
      (rr_max,method_max) = cbrt.calc_blended_rain_tropical(refl, zdr, kdp, rt_max, ice_flag=Ice_flag,
                                                            fhc=Fhc, predef=Predef, band=Band,
                                                            thresh_dz=Thresh_dz, thresh_zdr=Thresh_zdr,
                                                            thresh_kdp=Thresh_kdp,thresh_frac_ice=Thresh_frac_ice,
                                                            thresh_nexrad=Thresh_nexrad,
                                                            r_z_a=a_all, r_z_b=b_all,
                                                            r_z_a_c=a_conv, r_z_b_c=b_conv,
                                                            r_z_a_s=a_stra, r_z_b_s=b_stra,
                                                            r_kdp_a=a_kdp, r_kdp_b=b_kdp,
                                                            r_z_zdr_a=a_zh_zdr, r_z_zdr_b=b_ZH_zdr, r_z_zdr_c=b_zh_ZDR/10.,
                                                            fit_a=Fit_a, fit_b=Fit_b, method=Method,
                                                            r_kdp_zdr_a=a_kdp_zdr, r_kdp_zdr_b=b_KDP_zdr, r_kdp_zdr_c=b_kdp_ZDR/10.)

      # Add 2*RSME to all values
      rmse = np.zeros(rr_max.shape)
      for imethod in np.arange(1,len(zr_method)+1):
         if imethod == 6:
            mask = (method_max == imethod) & (rr_max < 10.)
            bin = rr_bin_stra['lt_10']
            rmse[mask] = rmse_a[imethod][bin] * (rr_max[mask] ** rmse_b[imethod][bin])

            mask = (method_max == imethod) & (rr_max > 20.)
            bin = rr_bin_stra['gt_20']
            rmse[mask] = rmse_a[imethod][bin] * (rr_max[mask] ** rmse_b[imethod][bin])
                    
            mask = (method_max == imethod) & (rr_max >= 10.) & (rr_max <= 20.) 
            bin = rr_bin_stra['10-20']
            rmse[mask] = rmse_a[imethod][bin] * (rr_max[mask] ** rmse_b[imethod][bin])
         else:
            mask = (method_max == imethod) & (rr_max < 20.)
            bin = rr_bin['lt_20']
            rmse[mask] = rmse_a[imethod][bin] * (rr_max[mask] ** rmse_b[imethod][bin])
                    
            mask = (method_max == imethod) & (rr_max > 60.)
            bin = rr_bin['gt_60']
            rmse[mask] = rmse_a[imethod][bin] * (rr_max[mask] ** rmse_b[imethod][bin])
                    
            mask = (method_max == imethod) & (rr_max >= 20.) & (rr_max <= 60.) 
            bin = rr_bin['20-60']
            rmse[mask] = rmse_a[imethod][bin] * (rr_max[mask] ** rmse_b[imethod][bin])
      rr_max = rr_max + (2*rmse)
                                                         
      # Add standard deviation of measurement error (sigma(epsilon_m)) to all values
      sigmaOverR = np.zeros(rr_max.shape)
      for imethod in np.arange(1,len(zr_method)+1):
         mask = (method_max == imethod)
         if imethod==3 or imethod==4 or imethod==5 or imethod==6:
            sigmaOverR[mask] = sigOverR[imethod]
         else:
            a = uncert_vals[imethod][uncert_vars['a']]
            b = uncert_vals[imethod][uncert_vars['b']]
            sigma_z = uncert_vals[imethod][uncert_vars['sigma_z']]
            sigma_zdr = uncert_vals[imethod][uncert_vars['sigma_zdr']]
            sigma_kdp = uncert_vals[imethod][uncert_vars['sigma_kdp']]
            if imethod==1:
               sigmaOverR[mask] = np.sqrt( (a**2) * (sigma_kdp**2)/(kdp[mask]**2) + 
                                           (b**2) * (sigma_zdr**2)/(zdr[mask]**2) )
            elif imethod==2:
               sigmaOverR[mask] = (a*sigma_kdp) * (a_kdp/rr_max[mask])**(1/a)
      rr_max = rr_max + (sigmaOverR * rr_max)
      rr_max[rr_max <= 0.0] = 0.0

                    
   # Reapply missing value masks to rr_min, rr_best and rr_max
   rr_min[rt_min_mask==1] = missing_value_rr
   rr_best[rt_best_mask==1] = missing_value_rr
   rr_max[rt_max_mask==1] = missing_value_rr
            
   return rr_min, rr_best, rr_max, method_min, method_best, method_max
