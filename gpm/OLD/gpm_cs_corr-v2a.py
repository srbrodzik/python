# Imports and functions
import os
import netCDF4 as nc4
import numpy as np
#import numpy.ma as ma

# User inputs
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/06/test'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/07/test'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/10/test'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/11/test'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/12/test'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2015/02/test'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2015/04/test'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2015/07/test'
indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2016/01/test'
THRESHOLD = 38.0  #dBZ
SLOPE_LO = 5.0
SLOPE_HI = 8.0
VERT_RES = 0.125  #km

for fname in os.listdir(indir):
    if fname.endswith('nc'):
        print fname

        ncid = nc4.Dataset(indir+'/'+fname,'a')

        refl = np.array(ncid.variables['refl'])
        (ntime,nalt,nlat,nlon) = refl.shape
        refl_missing = ncid.variables['refl'].missing_value
        bin_missing = int(refl_missing)

        rt = np.array(ncid.variables['rain_type'])
        rt_missing = ncid.variables['rain_type'].missing_value
        #rt[(rt==rt_missing)] = np.nan

        rt_uw = np.array(rt, copy=True)

        # change rain_type to stratiform in rt_uw if:
        # 1. rain_type = convective and
        # 2. max_refl in column >= THRESHOLD and
        # 3. (max_refl - dbz@(max_refl_ht-1km)/deltaHt > abs(SLOPE_LO)
        # 4. (max_refl - dbz@(max_refl_ht+1km)/deltaHt > abs(SLOPE_HI)

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

        # find bin & height of lowest good refl at each lat/lon
        min_ht_bin = np.zeros(rt.shape,dtype=int) + bin_missing
        for lev in range(0,nalt):
            min_ht_bin = np.where(np.logical_and( (refl[:,lev,:,:] != refl_missing),
                                                  (min_ht_bin == bin_missing) ), 
                                  lev, min_ht_bin ) 
        min_ht = min_ht_bin * (min_ht_bin != bin_missing) * VERT_RES
        min_ht[(min_ht_bin==bin_missing)] = refl_missing

        # Determine indices of interest 
        change = np.logical_and(rt==2,max_refl>=THRESHOLD)
        ind_time,ind_lat,ind_lon = np.where(change==True)
        npixels = ind_time.size

        # find slopes - check to make sure refl values are not missing
        for ind in range(0,npixels):
            itime = ind_time[ind]
            ilat = ind_lat[ind]
            ilon = ind_lon[ind]
            #print 'ind = ',ind,' and [itime,ilat,ilon] = [',itime,',',ilat,',',ilon,'] and max_refl = ',max_refl[itime,ilat,ilon]
            if max_refl[itime,ilat,ilon] != refl_missing:
                #print 'max_refl_bin = ',max_refl_bin[itime,ilat,ilon],' and max_refl = ',max_refl[itime,ilat,ilon]
                #print 'Looking for minus point'
                if max_refl_bin[itime,ilat,ilon] >= int(1.0/VERT_RES):
                    #print '   Not atypical'
                    ht_index = max_refl_bin[itime,ilat,ilon]-int(1.0/VERT_RES)
                    #print '   ht_index = ',ht_index,' and refl = ',refl[itime,ht_index,ilat,ilon]
                    refl_minus_diff = max_refl[itime,ilat,ilon]-refl[itime,ht_index,ilat,ilon]
                    ht_minus_diff = 1.0
                    #print '   refl_minus_diff = ',refl_minus_diff,' and ht_minus_diff = ',ht_minus_diff
                else:
                    #print '   Atypical'
                    ht_index = min_ht_bin[itime,ilat,ilon]
                    #print '   ht_index = ',ht_index,' and refl = ',refl[itime,ht_index,ilat,ilon]
                    refl_minus_diff = max_refl[itime,ilat,ilon]-refl[itime,ht_index,ilat,ilon]
                    ht_minus_diff = max_refl_ht[itime,ilat,ilon]-min_ht[itime,ilat,ilon]
                    #print '   refl_minus_diff = ',refl_minus_diff,' and ht_minus_diff = ',ht_minus_diff
                if ht_minus_diff != 0.0:
                    slope_minus = refl_minus_diff/ht_minus_diff
                else:
                    slope_minus = refl_missing
                #print '   slope_minus = ',slope_minus
    
                #print 'Looking for plus point'
                if max_ht_bin[itime,ilat,ilon] >= max_refl_bin[itime,ilat,ilon]+int(1.0/VERT_RES):
                    #print '   Not atypical'
                    ht_index = max_refl_bin[itime,ilat,ilon]+int(1.0/VERT_RES)
                    #print '   ht_index = ',ht_index,' and refl = ',refl[itime,ht_index,ilat,ilon]
                    refl_plus_diff = max_refl[itime,ilat,ilon]-refl[itime,ht_index,ilat,ilon]
                    ht_plus_diff = 1.0
                    #print '   refl_plus_diff = ',refl_plus_diff,' and ht_plus_diff = ',ht_plus_diff
                else:
                    #print '   Atypical'
                    ht_index = max_ht_bin[itime,ilat,ilon]
                    #print '   ht_index = ',ht_index,' and refl = ',refl[itime,ht_index,ilat,ilon]
                    refl_plus_diff = max_refl[itime,ilat,ilon]-refl[itime,ht_index,ilat,ilon]
                    #print 'refl = ',
                    ht_plus_diff = max_ht[itime,ilat,ilon]-max_refl_ht[itime,ilat,ilon]
                    #print '   refl_plus_diff = ',refl_plus_diff,' and ht_plus_diff = ',ht_plus_diff
                if ht_plus_diff != 0.0:
                    slope_plus = refl_plus_diff/ht_plus_diff
                else:
                    slope_plus = refl_missing
                #print '   slope_plus = ',slope_plus
            
                change = np.logical_and(abs(slope_minus) <= SLOPE_LO,
                                        abs(slope_plus) >= SLOPE_HI)
                #print 'change = ',change
                if change == True:
                    rt_uw[itime,ilat,ilon] = 1
                    
        # Add new raintype field to input netcdf file
        rt_uw_id = ncid.createVariable('rain_type_uw2','f4',
                                       ('time','latitude','longitude'),zlib=True)
        #rt_uw_id = ncid.createVariable('rain_type_uw','f4',
        #                              ('time','latitude','longitude'))
        rt_uw_id.units = "none"
        rt_uw_id.long_name = "Rain Type UW"
        rt_uw_id.stratiform = 1
        rt_uw_id.convective = 2
        rt_uw_id.other = 3
        rt_uw_id.missing_value = rt_missing
        rt_uw_id[:,:,:] = rt_uw
    
        ncid.close()

