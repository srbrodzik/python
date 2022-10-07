import os
import netCDF4 as nc4
import numpy as np
#import numpy.ma as ma

# User inputs
indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/06/test'
THRESHOLD = 40.0  #dBZ
MIN_HT = 3.0      #km
MAX_HT = 5.0      #km
VERT_RES = 0.125  #km
#fltr_missing = -99

for fname in os.listdir(indir):

    if fname.endswith('nc'):

        print fname

        # open fname
        ncid = nc4.Dataset(indir+'/'+fname,'a')
        
        # read vars of interest (refl and rain_type)
        refl = ncid.variables['refl'][:]
        (ntime,nalt,nlat,nlon) = refl.shape
        rt = ncid.variables['rain_type'][:]
        rt_missing = ncid.variables['rain_type'].missing_value
        srt = ncid.variables['shallow_rain_type'][:]
        #print 'Done reading vars of interest'
        
        # create arrays for rain_type_uw and filters
        rt_uw = np.array(rt, copy=True)
        #fltr = ma.masked_array(np.zeros(rt.shape), rt.mask)
        #ma.set_fill_value(fltr, fltr_missing)
        #fltr_srt = ma.masked_array(np.zeros(rt.shape), rt.mask)
        #ma.set_fill_value(fltr_srt, fltr_missing)
        #fltr_total = ma.masked_array(np.zeros(rt.shape), rt.mask)
        #ma.set_fill_value(fltr_total, fltr_missing)
        #print 'Done creating output arrays'
    
        # change rain_type values to STRATIFORM in rt_uw if:
        # 1. rain_type = CONVECTIVE and
        # 2. max refl in column > THRESHOLD and
        # 3. ht of max refl in column is between MIN_HT and MAX_HT
        conv_mask = np.zeros(rt.shape,dtype=int)
        conv_mask[(rt==2)] = 1
        max_refl = np.zeros(rt.shape,dtype=float)
        max_refl = np.amax(refl,axis=1)
        max_refl_ht = np.zeros(rt.shape,dtype=float)
        max_refl_ht = np.argmax(refl,axis=1)

        # find gridpoints where rain_type will be changed from CONV to STRA
        chg = np.logical_and(np.logical_and(conv_mask==1,max_refl>=THRESHOLD),
                             np.logical_and(max_refl_ht*VERT_RES>=MIN_HT,
                                            max_refl_ht*VERT_RES<=MAX_HT))
        rt_uw[chg==True]=1
        #fltr[chg == True] = 1
        #print 'Done changing rain_type values'

        # find gridpoints where rain is shallow (isolated or non-isolated)
        # with high certainty
        #  0 => no_shallow_rain
        # 10 => shallow_isolated_maybe
        # 11 => shallow_isolated_certain
        # 20 => shallow_nonisolated_maybe
        # 21 => shallow_nonisolated_certain
        #chg_srt = np.logical_or(srt==11,srt==21)
        #fltr_srt[chg_srt == True] = 1
        #print 'Done creating srt mask'

        # find gridpoints where rain_type will change and where rain is shallow
        #chg_total = np.logical_and(chg,chg_srt)
        #fltr_total[chg_total == True] = 1
        #print 'Done creating total mask'

        # add new fields to netcdf file
        rt_uw_id = ncid.createVariable('rain_type_uw','f4',
                                       ('time','latitude','longitude'),zlib=True)
        rt_uw_id.units = "none"
        rt_uw_id.long_name = "Rain Type UW"
        rt_uw_id.stratiform = 1
        rt_uw_id.convective = 2
        rt_uw_id.other = 3
        rt_uw_id.missing_value = rt_missing
        rt_uw_id[:,:,:] = rt_uw

        #filter_uw_id = ncid.createVariable('filter_chg_rt','f4',
        #                                   ('time','latitude','longitude'),zlib=True)
        #filter_uw_id.units = "none"
        #filter_uw_id.long_name = "Rain type change mask"
        #filter_uw_id.missing_value = fltr_missing
        #filter_uw_id[:,:,:] = fltr
    
        #filter_srt_id = ncid.createVariable('filter_srt','f4',
        #                                    ('time','latitude','longitude'),zlib=True)
        #filter_srt_id.units = "none"
        #filter_srt_id.long_name = "Shallow rain mask"
        #filter_srt_id.missing_value = fltr_missing
        #filter_srt_id[:,:,:] = fltr_srt
    
        #filter_tot_id = ncid.createVariable('filter_total','f4',
        #                                    ('time','latitude','longitude'),zlib=True)
        #filter_tot_id.units = "none"
        #filter_tot_id.long_name = "Rain type change and shallow rain mask"
        #filter_tot_id.missing_value = fltr_missing
        #filter_tot_id[:,:,:] = fltr_total
        #print 'Done adding fields to netcdf file'
        
        # close input file
        ncid.close()
        #print 'Done closing netcdf file'


