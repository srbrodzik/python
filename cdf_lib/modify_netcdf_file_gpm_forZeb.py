import os
import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

#NOTE: Use ~/python/gpm/findWCCfiles.py to copy files to indir

# gpm_ku_v05_uw_ashwin
#indir = '/home/disk/radar/india_main/zeb-india-data/gpm_ku_ashwin/netcdf4'
#outdir = '/home/disk/radar/india_main/zeb-india-data/gpm_ku_ashwin/netcdfZeb'

# gpm_ku_v05_uw
#indir  = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str'
#outdir = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_wcc_str'
#indir  = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dwc_str'
#outdir = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_dwc_str'

# gpm_ku_ite
#indir  = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_ite/netcdf4_wcc_str'
#outdir = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_ite/netcdfZeb_wcc_str'
#indir  = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_ite/netcdf4_dwc_str'
#outdir = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_ite/netcdfZeb_dwc_str'

# gpm_ku_asia
#indir  = '/home/disk/bob/gpm/asia_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dcc_str'
#outdir = '/home/disk/bob/gpm/asia_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_dcc_str'
#indir  = '/home/disk/bob/gpm/asia_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str'
#outdir = '/home/disk/bob/gpm/asia_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_wcc_str'

# gpm_ku_nam
#indir  = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_bsr_str_mask'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_bsr_str_mask'
#indir  = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str_mask'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_wcc_str_mask'

# gpm_ku_aka
#indir  = '/home/disk/bob/gpm/aka_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dcc_str'
#outdir = '/home/disk/bob/gpm/aka_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_dcc_str'

# gpm_ku_eur
#indir  = '/home/disk/bob/gpm/eur_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dcc_str'
#outdir = '/home/disk/bob/gpm/eur_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_dcc_str'

# gpm_ku_h02
#indir  = '/home/disk/bob/gpm/h02_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_mod'
#outdir = '/home/disk/bob/gpm/h02_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_wcc_mod'

# gpm_ku_cio
#indir  = '/home/disk/bob/gpm/cio_ku/zeb-data/gpm_ku_v05_uw/netcdf4'
#outdir = '/home/disk/bob/gpm/cio_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb'

# gpm_ku_epo
#indir  = '/home/disk/bob/gpm/epo_ku/zeb-data/gpm_ku_v05_uw/netcdf4'
#outdir = '/home/disk/bob/gpm/epo_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb'

# gpm_ku_sam_rel
indir  = '/home/disk/bob/gpm/sam_rel_ku/classify/ex_data_v06/2019/04'
outdir = '/home/disk/monsoon/relampago/zeb/gpm_ku'

#uw = 'true'
uw = 'false'

masks = 'true'
#masks = 'false'

thres = 'str'
#thres = 'mod'
#thres = 'both'

for fname in os.listdir(indir):
    print fname
    if fname.endswith('nc'):
        
        # open compressed netcdf4 file for reading
        ncid = nc4.Dataset(indir+'/'+fname,'r')

        # read vars
        times = ncid.variables['time'][:]
        ntimes = times.size
        lats = ncid.variables['lat'][:]
        nlats = lats.size
        lons = ncid.variables['lon'][:]
        nlons = lons.size
        alts = ncid.variables['alt'][:]
        nalts = alts.size
        nalts_125 = 125
        ntimes = 1
        time = ncid.variables['time'][:]
        refl = ncid.variables['refl'][:]
        refl_missing = ncid.variables['refl']._FillValue
        rtype = ncid.variables['rain_type'][:]
        rtype_missing = ncid.variables['rain_type']._FillValue
        if uw == 'true':
            rtype_uw = ncid.variables['rain_type_uw'][:]
            rtype_uw_missing = ncid.variables['rain_type_uw']._FillValue
        #rtype_uw2 = ncid.variables['rain_type_uw2'][:]
        #rtype_uw2_missing = ncid.variables['rain_type_uw2']._FillValue
        srt = ncid.variables['shallow_rain_type'][:]
        srt_missing = ncid.variables['shallow_rain_type']._FillValue
        precip = ncid.variables['near_surf_rain'][:]
        precip_missing = ncid.variables['near_surf_rain']._FillValue
        swath = ncid.variables['swath'][:]
        swath_missing = ncid.variables['swath']._FillValue
        if masks == 'true':
            if thres == 'str' or thres == 'both':
                bsr_str = ncid.variables['bsr_mask_str'][:]
                bsr_str_missing = ncid.variables['bsr_mask_str']._FillValue
                dcc_str = ncid.variables['dcc_mask_str'][:]
                dcc_str_missing = ncid.variables['dcc_mask_str']._FillValue
                dwc_str = ncid.variables['dwc_mask_str'][:]
                dwc_str_missing = ncid.variables['dwc_mask_str']._FillValue
                wcc_str = ncid.variables['wcc_mask_str'][:]
                wcc_str_missing = ncid.variables['wcc_mask_str']._FillValue
                storm_str = ncid.variables['storm_mask_str'][:]
                storm_str_missing = ncid.variables['storm_mask_str']._FillValue
                #replace all positive storm values with 1
                storm_str[storm_str>0] = 1
            if thres == 'mod' or thres == 'both':
                bsr_mod = ncid.variables['bsr_mask_mod'][:]
                bsr_mod_missing = ncid.variables['bsr_mask_mod']._FillValue
                dcc_mod = ncid.variables['dcc_mask_mod'][:]
                dcc_mod_missing = ncid.variables['dcc_mask_mod']._FillValue
                dwc_mod = ncid.variables['dwc_mask_mod'][:]
                dwc_mod_missing = ncid.variables['dwc_mask_mod']._FillValue
                wcc_mod = ncid.variables['wcc_mask_mod'][:]
                wcc_mod_missing = ncid.variables['wcc_mask_mod']._FillValue
                storm_mod = ncid.variables['storm_mask_mod'][:]
                storm_mod_missing = ncid.variables['storm_mask_mod']._FillValue
                #replace all positive storm values with 1
                storm_mod[storm_mod>0] = 1

        # read global attributes
        title = ncid.title
        orbit = ncid.orbit
        lat_min = ncid.lat_min
        lat_max = ncid.lat_max
        lon_min = ncid.lon_min
        lon_max = ncid.lon_max

        # close input file
        ncid.close()

        # create new uncompressed netcdf file with subset of fields for zeb analysis
        ##parts = fname.split('_')
        ##ncname = outdir+'/'+parts[0]+'_'+parts[1]+'_'+parts[2]+'_'+parts[5]+'_'+parts[6]
        ncname = outdir+'/'+fname
        # NOTES:
        #   NETCDF3_CLASSIC - limited to file sizes less than 2GB
        #   NETCDF3_64BIT_OFFSET - allows for file sizes greater than 2GB
        #   NETCDF3_64BIT_DATA - extends NETCDF3_64BIT_OFFSET to allow for unsigned/64 bit int data types & 64-bit dim sizes
        #   NETCDF4 - default
        ncid = nc4.Dataset(ncname,'w',format='NETCDF3_CLASSIC')

        # create dims
        time = ncid.createDimension('time',None) # None implies UNLIMITED
        z0 = ncid.createDimension('altitude',nalts_125)
        y0 = ncid.createDimension('latitude',nlats)
        x0 = ncid.createDimension('longitude',nlons)

        # create vars
        base_time_var = ncid.createVariable('base_time',np.int32)
        time_offset_var = ncid.createVariable('time_offset',np.float32,('time') )
        #jday_var = ncid.createVariable('julian_date',np.float64,('time') )
        z0_var = ncid.createVariable('altitude',np.float32,('altitude') )
        y0_var = ncid.createVariable('latitude',np.float32,('latitude') )
        x0_var = ncid.createVariable('longitude',np.float32,('longitude') )
        refl_var = ncid.createVariable('refl',np.float32,('time','altitude','latitude','longitude'),
                                       fill_value=refl_missing )
        rtype_var = ncid.createVariable('rain_type',np.float32,('time','latitude','longitude'),
                                        fill_value=rtype_missing )
        if uw == 'true':
            rtype_uw_var = ncid.createVariable('rain_type_uw',np.float32,('time','latitude','longitude'),
                                               fill_value=rtype_uw_missing )
        srt_var = ncid.createVariable('shallow_rain_type',np.float32,('time','latitude','longitude'),
                                      fill_value=srt_missing )
        precip_var = ncid.createVariable('near_surf_rain',np.float32,('time','latitude','longitude'),
                                         fill_value=precip_missing )
        swath_var = ncid.createVariable('swath',np.float32,('time','latitude','longitude'),
                                        fill_value=swath_missing )
        if masks == 'true':
            if thres == 'str' or thres == 'both':
                bsr_str_var = ncid.createVariable('bsr_mask_str',np.float32,('time','latitude','longitude'),
                                              fill_value=bsr_str_missing )
                dcc_str_var = ncid.createVariable('dcc_mask_str',np.float32,('time','latitude','longitude'),
                                              fill_value=dcc_str_missing )
                dwc_str_var = ncid.createVariable('dwc_mask_str',np.float32,('time','latitude','longitude'),
                                              fill_value=dwc_str_missing )
                wcc_str_var = ncid.createVariable('wcc_mask_str',np.float32,('time','latitude','longitude'),
                                              fill_value=wcc_str_missing )
                storm_str_var = ncid.createVariable('storm_mask_str',np.float32,('time','latitude','longitude'),
                                                fill_value=storm_str_missing )
            if thres == 'mod' or thres == 'both':
                bsr_mod_var = ncid.createVariable('bsr_mask_mod',np.float32,('time','latitude','longitude'),
                                              fill_value=bsr_mod_missing )
                dcc_mod_var = ncid.createVariable('dcc_mask_mod',np.float32,('time','latitude','longitude'),
                                              fill_value=dcc_mod_missing )
                dwc_mod_var = ncid.createVariable('dwc_mask_mod',np.float32,('time','latitude','longitude'),
                                              fill_value=dwc_mod_missing )
                wcc_mod_var = ncid.createVariable('wcc_mask_mod',np.float32,('time','latitude','longitude'),
                                              fill_value=wcc_mod_missing )
                storm_mod_var = ncid.createVariable('storm_mask_mod',np.float32,('time','latitude','longitude'),
                                                fill_value=storm_mod_missing )

        # define attributes
        base_time_var.standard_name = 'time'
        base_time_var.units = 'seconds since 1970-01-01 00:00:00 +0000'

        time_offset_var.standard_name = 'time offset'
        time_offset_var.units = 'seconds since base_time'

        #jday_var.long_name = 'Julian Date Number'
        #jday_var.units = 'Julian number after October 4, 1582-10-4 12:00:00'

        z0_var.units = 'km'
        #z0_var.positive = 'up'
        #z0_var.axis = 'Z'

        #y0_var.standard_name = 'latitude'
        y0_var.long_name = 'Latitude coordinate'
        y0_var.units = 'degrees_north'
        #y0_var.axis = 'Y'

        #x0_var.standard_name = 'longitude'
        x0_var.long_name = 'Longitude coordinate'
        x0_var.units = 'degrees_west'
        #x0_var.axis = 'X'

        refl_var.long_name = 'GPM-Ku Reflectivity'
        refl_var.units = 'dBZ'
        refl_var.missing_value = refl_missing

        rtype_var.long_name = 'Rain Type'
        rtype_var.units = 'none'
        rtype_var.stratiform = 1
        rtype_var.convective = 2
        rtype_var.other = 3
        rtype_var.missing_value = rtype_missing

        if uw == 'true':
            rtype_uw_var.long_name = 'Rain Type UW'
            rtype_uw_var.units = 'none'
            rtype_uw_var.stratiform = 1
            rtype_uw_var.convective = 2
            rtype_uw_var.other = 3
            rtype_uw_var.missing_value = rtype_uw_missing

        srt_var.long_name = 'Shallow Rain Flag'
        srt_var.units = 'none'
        srt_var.no_shallow_rain = 0
        srt_var.shallow_isolated_maybe = 10
        srt_var.shallow_isolated_certain = 11
        srt_var.shallow_nonisolated_maybe = 20
        srt_var.shallow_nonisolated_certain = 21
        srt_var.no_rain_value = -1111
        srt_var.missing_value = srt_missing
        
        precip_var.long_name = 'Near Surface Rain'
        precip_var.units = 'mm/hr'
        precip_var.missing_value = precip_missing

        swath_var.long_name = 'GPM-Ku coverage area'
        swath_var.units = 'none'
        swath_var.swath = 0
        swath_var.no_swath = 1
        swath_var.missing_value = swath_missing

        if masks == 'true':
            if thres == 'str' or thres == 'both':
                bsr_str_var.long_name = 'BSR core mask strong'
                bsr_str_var.units = 'none'
                bsr_str_var.missing_value = bsr_str_missing

                dcc_str_var.long_name = 'DCC core mask strong'
                dcc_str_var.units = 'none'
                dcc_str_var.missing_value = dcc_str_missing

                dwc_str_var.long_name = 'DWC core mask strong'
                dwc_str_var.units = 'none'
                dwc_str_var.missing_value = dwc_str_missing

                wcc_str_var.long_name = 'WCC core mask strong'
                wcc_str_var.units = 'none'
                wcc_str_var.missing_value = wcc_str_missing

                storm_str_var.long_name = 'Storm mask strong'
                storm_str_var.units = 'none'
                storm_str_var.missing_value = storm_str_missing
            if thres == 'mod' or thres == 'both':
                bsr_mod_var.long_name = 'BSR core mask moderate'
                bsr_mod_var.units = 'none'
                bsr_mod_var.missing_value = bsr_mod_missing

                dcc_mod_var.long_name = 'DCC core mask moderate'
                dcc_mod_var.units = 'none'
                dcc_mod_var.missing_value = dcc_mod_missing

                dwc_mod_var.long_name = 'DWC core mask moderate'
                dwc_mod_var.units = 'none'
                dwc_mod_var.missing_value = dwc_mod_missing

                wcc_mod_var.long_name = 'WCC core mask moderate'
                wcc_mod_var.units = 'none'
                wcc_mod_var.missing_value = wcc_mod_missing

                storm_mod_var.long_name = 'Storm mask moderate'
                storm_mod_var.units = 'none'
                storm_mod_var.missing_value = storm_mod_missing

        # define global attributes
        ncid.title = title
        ncid.orbit = orbit
        ncid.lat_min = lat_min
        ncid.lat_max = lat_max
        ncid.lon_min = lon_min
        ncid.lon_max = lon_max

        # write vars to file
        base_time_var[:] = times
        time_offset_var[:] = [0]
        #jday_var[:] = jdate
        z0_var[:] = alts[0:125]
        y0_var[:] = lats
        x0_var[:] = lons
        refl_var[:,:,:,:] = refl[:,0:125,:,:]
        rtype_var[:,:,:] = rtype
        if uw == 'true':
            rtype_uw_var[:,:,:] = rtype_uw
        srt_var[:,:,:] = srt
        precip_var[:,:,:] = precip
        swath_var[:,:,:] = swath
        if masks == 'true':
            if thres == 'str' or thres == 'both':
                bsr_str_var[:,:,:] = bsr_str
                dcc_str_var[:,:,:] = dcc_str
                dwc_str_var[:,:,:] = dwc_str
                wcc_str_var[:,:,:] = wcc_str
                storm_str_var[:,:,:] = storm_str
            if thres == 'mod' or thres == 'both':
                bsr_mod_var[:,:,:] = bsr_mod
                dcc_mod_var[:,:,:] = dcc_mod
                dwc_mod_var[:,:,:] = dwc_mod
                wcc_mod_var[:,:,:] = wcc_mod
                storm_mod_var[:,:,:] = storm_mod

        # close file
        ncid.close()


