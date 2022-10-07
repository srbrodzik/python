# Imports and functions
import os
import netCDF4 as nc4
import numpy as np
#import numpy.ma as ma

# For plotting
#%matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

# User inputs
indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05_testing'
years = ['2016']
months = ['12']
VERT_RES = 0.125  #km
STRA = 1
CONV = 2

# ORIG
#THRESHOLD = 38.0  #dBZ
#DELTA_HT_FOR_SLOPE = 1.5  #km
#SLOPE = 8.0
#MAX_BB_HT = 5.0   #km
#REFL_MAX_OFFSET = 2.0 #dBZ

# TRY1
#THRESHOLD = 0.0  #dBZ
#DELTA_HT_FOR_SLOPE = 1.0  #km
#SLOPE = 5.0
#MAX_BB_HT = 5.0   #km
#REFL_MAX_OFFSET = 2.0 #dBZ

# TRY2
THRESHOLD = 0.0  #dBZ
DELTA_HT_FOR_SLOPE = 1.5  #km
SLOPE = 7.5
MAX_BB_HT = 5.0   #km
REFL_MAX_OFFSET = 2.0 #dBZ

#

for iyear in years:
    for imonth in months:
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
            if fname.endswith('nc'):
                print fname

                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'a')

                refl = np.array(ncid.variables['refl'])
                refl = refl[:,1:,:,:]  # get rid of lowest level due to suspicious data
                (ntime,nalt,nlat,nlon) = refl.shape
                #refl_missing = ncid.variables['refl'].missing_value
                refl_missing = ncid.variables['refl']._FillValue
                bin_missing = int(refl_missing)

                rt = np.array(ncid.variables['rain_type'])
                #rt_missing = ncid.variables['rain_type'].missing_value
                rt_missing = ncid.variables['rain_type']._FillValue
                #rt[(rt==rt_missing)] = np.nan

                rt_uw = np.array(rt, copy=True)

                srt = np.array(ncid.variables['shallow_rain_type'])
                #srt_missing = ncid.variables['shallow_rain_type'].missing_value
                srt_missing = ncid.variables['shallow_rain_type']._FillValue
                #srt[(srt==srt_missing)] = np.nan

                # find max_refl and max_refl_ht at each lat/lon
                max_refl = np.zeros(rt.shape,dtype=float)+refl_missing
                max_refl = np.amax(refl,axis=1)
                max_refl_bin = np.zeros(rt.shape,dtype=int)+bin_missing
                max_refl_bin = np.argmax(refl,axis=1)
                max_refl_ht = max_refl_bin * VERT_RES

                max_refl_bin[(max_refl==refl_missing)] = bin_missing
                max_refl_ht[(max_refl==refl_missing)] = refl_missing

                # initialize reference refl/bin/ht
                # these reflect highest refl,bin,ht above max_refl,bin.ht and within +/- REFL_MAX_OFFSET
                reference_refl = np.array(max_refl, copy=True)
                reference_refl_bin = np.array(max_refl_bin, copy=True)
                reference_refl_ht = np.array(max_refl_ht, copy=True)
                
                # find bin & height of highest good refl at each lat/lon
                max_ht_bin = np.zeros(rt.shape,dtype=int) + bin_missing
                for lev in range(nalt-1,-1,-1):
                    max_ht_bin = np.where(np.logical_and((refl[:,lev,:,:] != refl_missing),
                                                         (max_ht_bin == bin_missing) ), 
                                          lev, max_ht_bin )
                max_ht = max_ht_bin * (max_ht_bin != bin_missing) * VERT_RES
                max_ht[(max_ht_bin==bin_missing)] = refl_missing

                # For rt_uw, change rain_type to stratiform if:
                # 1. rain_type = convective and
                # 2. shallow_rain_type = 0 and
                # 3. max_refl in column >= THRESHOLD and   ## DON'T USE THIS CRITERIA ANYMORE
                # 4. [(dbz @ reference point) - (dbz @ point 1.5km higher)]/
                #    [(ht  @ reference point) - (ht  @ point 1.5km higher)] >= abs(SLOPE)
                #  as long as ht @ reference point <= MAX_BB_HT
                
                # Determine indices of interest 
                change = np.logical_and(np.logical_and(rt==2,srt==0),max_refl>=THRESHOLD)
                ind_time,ind_lat,ind_lon = np.where(change==True)
                npixels = ind_time.size

                # for each index
                #    start at max_refl_ht & go up to find reference refl within +/- REFL_MAX_OFFSET
                #    if height of new reference refl is less than MAX_BB_HT
                #        determine slope between that point and point DELTA_HT_FOR_SLOPE higher or highest good point in column
                #        if slope exceeds SLOPE, change raintype from convective to stratiform                    
                for ind in range(0,npixels):
                    itime = ind_time[ind]
                    ilat = ind_lat[ind]
                    ilon = ind_lon[ind]
                    
                    if max_refl[itime,ilat,ilon] != refl_missing:
                        print('Max refl point: ',max_refl[itime,ilat,ilon],max_refl_bin[itime,ilat,ilon],max_refl_ht[itime,ilat,ilon])
                        
                        next_bin = reference_refl_bin[itime,ilat,ilon]+1
                        next_ht = reference_refl_ht[itime,ilat,ilon]+VERT_RES
                        next_refl = refl[itime,next_bin,ilat,ilon]
                        #print(reference_refl[itime,ilat,ilon],reference_refl_bin[itime,ilat,ilon],reference_refl_ht[itime,ilat,ilon])
                        #print(next_refl,next_bin,next_ht)
                        while abs(max_refl[itime,ilat,ilon]-next_refl) <= REFL_MAX_OFFSET:
                            reference_refl[itime,ilat,ilon]=next_refl
                            reference_refl_bin[itime,ilat,ilon] = next_bin
                            reference_refl_ht[itime,ilat,ilon]=next_ht  
                            next_bin = reference_refl_bin[itime,ilat,ilon]+1
                            next_ht = reference_refl_ht[itime,ilat,ilon]+VERT_RES
                            next_refl = refl[itime,next_bin,ilat,ilon]
                            #print(reference_refl[itime,ilat,ilon],reference_refl_bin[itime,ilat,ilon],reference_refl_ht[itime,ilat,ilon])
                            #print(next_refl,next_bin,next_ht)

                        print('Reference point: ',reference_refl[itime,ilat,ilon],reference_refl_bin[itime,ilat,ilon],reference_refl_ht[itime,ilat,ilon])
                        
                        if reference_refl_ht[itime,ilat,ilon] < MAX_BB_HT:
                            if (reference_refl_ht[itime,ilat,ilon] + DELTA_HT_FOR_SLOPE) <= max_ht[itime,ilat,ilon]:
                                compare_bin = reference_refl_bin[itime,ilat,ilon] + int(DELTA_HT_FOR_SLOPE/VERT_RES)
                                compare_ht = reference_refl_ht[itime,ilat,ilon] + DELTA_HT_FOR_SLOPE
                                compare_refl = refl[itime,compare_bin,ilat,ilon]
                            else:
                                compare_bin = max_ht_bin[itime,ilat,ilon]
                                compare_ht = max_ht[itime,ilat,ilon]
                                compare_refl = refl[itime,compare_bin,ilat,ilon]
                            print('Compare point: ',compare_refl,compare_bin,compare_ht)
                            if ( abs(reference_refl[itime,ilat,ilon] - compare_refl)/
                                 abs(reference_refl_ht[itime,ilat,ilon] - compare_ht) ) >= SLOPE:
                                rt_uw[itime,ilat,ilon] = STRA

                # Add new raintype field to input netcdf file
                rt_uw_id = ncid.createVariable('rain_type_uw2','f4',
                                               #('time','lat','lon'),zlib=True)
                                               ('time','lat','lon'),zlib=True,fill_value=rt_missing)
                rt_uw_id.units = "none"
                rt_uw_id.long_name = "Rain Type UW 2"
                rt_uw_id.stratiform = 1
                rt_uw_id.convective = 2
                rt_uw_id.other = 3
                #rt_uw_id.missing_value = rt_missing
                rt_uw_id[:,:,:] = rt_uw
    
                ncid.close()

