# Imports and functions
import os
import netCDF4 as nc4
import numpy as np
#import numpy.ma as ma

# User inputs
indir = '/home/disk/bob/gpm/asia_ku/classify/ex_data_v05'
#years = ['2014']
#months = ['03','04','05','06','07','08','09','10','11','12']
#years = ['2015','2016']
#months = ['01','02','03','04','05','06','07','08','09','10','11','12']
years = ['2017']
months = ['01','02','03','04','05','06','07','08']
THRESHOLD = 38.0  #dBZ
SLOPE_HI = 8.0
VERT_RES = 0.125  #km

for iyear in years:
    for imonth in months:
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
            if fname.endswith('nc'):
                print fname

                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'a')

                refl = np.array(ncid.variables['refl'])
                (ntime,nalt,nlat,nlon) = refl.shape
                #refl_missing = ncid.variables['refl'].missing_value
                refl_missing = ncid.variables['refl']._FillValue
                bin_missing = int(refl_missing)

                rt = np.array(ncid.variables['rain_type'])
                #rt_missing = ncid.variables['rain_type'].missing_value
                rt_missing = ncid.variables['rain_type']._FillValue
                #rt[(rt==rt_missing)] = np.nan

                rt_uw = np.array(rt, copy=True)

                # For rt_uw, change rain_type to stratiform if:
                # 1. rain_type = convective and
                # 2. max_refl in column >= THRESHOLD and
                # 3. (max_refl - dbz@(max_refl_ht+1km)/deltaHt >= abs(SLOPE_HI)    

                # find max_refl and max_refl_ht at each lat/lon
                max_refl = np.zeros(rt.shape,dtype=float)+refl_missing
                max_refl = np.amax(refl,axis=1)
                max_refl_bin = np.zeros(rt.shape,dtype=int)+bin_missing
                max_refl_bin = np.argmax(refl,axis=1)
                max_refl_ht = max_refl_bin * VERT_RES

                max_refl_bin[(max_refl==refl_missing)] = bin_missing
                max_refl_ht[(max_refl==refl_missing)] = refl_missing

                # find bin & height of highest good refl at each lat/lon
                max_ht_bin = np.zeros(rt.shape,dtype=int) + bin_missing
                for lev in range(nalt-1,-1,-1):
                    max_ht_bin = np.where(np.logical_and((refl[:,lev,:,:] != refl_missing),
                                                         (max_ht_bin == bin_missing) ), 
                                          lev, max_ht_bin )
                max_ht = max_ht_bin * (max_ht_bin != bin_missing) * VERT_RES
                max_ht[(max_ht_bin==bin_missing)] = refl_missing

                # Determine indices of interest 
                change = np.logical_and(rt==2,max_refl>=THRESHOLD)
                ind_time,ind_lat,ind_lon = np.where(change==True)
                npixels = ind_time.size

                # find slope hi - check to make sure refl values are not missing
                for ind in range(0,npixels):
                    itime = ind_time[ind]
                    ilat = ind_lat[ind]
                    ilon = ind_lon[ind]
                    if max_refl[itime,ilat,ilon] != refl_missing:
                        if max_ht_bin[itime,ilat,ilon] >= max_refl_bin[itime,ilat,ilon]+int(1.5/VERT_RES):
                            ht_index = max_refl_bin[itime,ilat,ilon]+int(1.5/VERT_RES)
                            refl_plus_diff = max_refl[itime,ilat,ilon]-refl[itime,ht_index,ilat,ilon]
                            ht_plus_diff = 1.5
                        else:
                            ht_index = max_ht_bin[itime,ilat,ilon]
                            refl_plus_diff = max_refl[itime,ilat,ilon]-refl[itime,ht_index,ilat,ilon]
                            ht_plus_diff = max_ht[itime,ilat,ilon]-max_refl_ht[itime,ilat,ilon]
                        if ht_plus_diff != 0.0:
                            slope_plus = refl_plus_diff/ht_plus_diff
                        else:
                            slope_plus = refl_missing
            
                        change3 = abs(slope_plus) >= SLOPE_HI
                        #print '   change3 = ',change3
                        if change3 == True:
                            rt_uw[itime,ilat,ilon] = 1
                    
                # Add new raintype field to input netcdf file
                rt_uw_id = ncid.createVariable('rain_type_uw','f4',
                                               #('time','lat','lon'),zlib=True)
                                               ('time','lat','lon'),zlib=True,fill_value=rt_missing)
                rt_uw_id.units = "none"
                rt_uw_id.long_name = "Rain Type UW"
                rt_uw_id.stratiform = 1
                rt_uw_id.convective = 2
                rt_uw_id.other = 3
                #rt_uw_id.missing_value = rt_missing
                rt_uw_id[:,:,:] = rt_uw
    
                ncid.close()

