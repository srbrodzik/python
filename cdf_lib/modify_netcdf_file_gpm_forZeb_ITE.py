import os
import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

#NOTE: Use ~/python/gpm/findWCCfiles.py to copy files to indir

#indir  = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_ite/netcdf4_xxx_str'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_ite/netcdfZeb_xxx_str'

#indir  = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dwc_str'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdfZeb_dwc_str'
indir  = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_ite/netcdf4_dwc_str'
outdir = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_ite/netcdfZeb_dwc_str'

uw = 'false'
useCores = 'false'

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
        if useCores == 'true':
            bsr = ncid.variables['bsr_mask_str'][:]
            bsr_missing = ncid.variables['bsr_mask_str']._FillValue
            dcc = ncid.variables['dcc_mask_str'][:]
            dcc_missing = ncid.variables['dcc_mask_str']._FillValue
            dwc = ncid.variables['dwc_mask_str'][:]
            dwc_missing = ncid.variables['dwc_mask_str']._FillValue
            wcc = ncid.variables['wcc_mask_str'][:]
            wcc_missing = ncid.variables['wcc_mask_str']._FillValue
            storm = ncid.variables['storm_mask_str'][:]
            storm_missing = ncid.variables['storm_mask_str']._FillValue
            #replace all positive storm values with 1
            storm[storm>0] = 1

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
        if useCores == 'true':
            bsr_var = ncid.createVariable('bsr_mask_str',np.float32,('time','latitude','longitude'),
                                          fill_value=bsr_missing )
            dcc_var = ncid.createVariable('dcc_mask_str',np.float32,('time','latitude','longitude'),
                                          fill_value=dcc_missing )
            dwc_var = ncid.createVariable('dwc_mask_str',np.float32,('time','latitude','longitude'),
                                          fill_value=dwc_missing )
            wcc_var = ncid.createVariable('wcc_mask_str',np.float32,('time','latitude','longitude'),
                                          fill_value=wcc_missing )
            storm_var = ncid.createVariable('storm_mask_str',np.float32,('time','latitude','longitude'),
                                            fill_value=storm_missing )

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

        if useCores == 'true':
        
            bsr_var.long_name = 'BSR core mask strong'
            bsr_var.units = 'none'
            bsr_var.missing_value = bsr_missing

            dcc_var.long_name = 'DCC core mask strong'
            dcc_var.units = 'none'
            dcc_var.missing_value = dcc_missing

            dwc_var.long_name = 'DWC core mask strong'
            dwc_var.units = 'none'
            dwc_var.missing_value = dwc_missing

            wcc_var.long_name = 'WCC core mask strong'
            wcc_var.units = 'none'
            wcc_var.missing_value = wcc_missing

            storm_var.long_name = 'Storm mask strong'
            storm_var.units = 'none'
            storm_var.missing_value = storm_missing

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
        if useCores == 'true':
            bsr_var[:,:,:] = bsr
            dcc_var[:,:,:] = dcc
            dwc_var[:,:,:] = dwc
            wcc_var[:,:,:] = wcc
            storm_var[:,:,:] = storm

        # close file
        ncid.close()


