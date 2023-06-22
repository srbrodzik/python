import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt
import sys

def writeCFnetcdf(ncname,core_thresh,method,a,b,tune_thresh,sm_rad,
                  fill_dbz,bg_diff,institution,source,title,
                  references,comment,time_val,lat_vals,lon_vals,lat_origin,
                  lon_origin,rain_type_vals,rain_type_basic_vals,
                  conv_core_vals,bkgrnd_vals,missing_value):
    
    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert time to date and time
    date = dt.datetime.utcfromtimestamp(time_val[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    ncid.createDimension('time',None) # None implies UNLIMITED
    ncid.createDimension('lon',lon_vals.shape[0])
    ncid.createDimension('lat',lat_vals.shape[0])
    ncid.createDimension('nchar', 4)

    # create variables
    timeVar = ncid.createVariable('time',np.float64,('time'),zlib=True )
    latVar = ncid.createVariable('lat',np.float32,('lat'),zlib=True )
    lonVar = ncid.createVariable('lon',np.float32,('lon'),zlib=True )
    #latOrigVar = ncid.createVariable('lat_origin',np.float32,zlib=True )
    #lonOrigVar = ncid.createVariable('lon_origin',np.float32,zlib=True )
    ctVar = ncid.createVariable('core_thresh',np.float32)
    #methVar = ncid.createVariable('method',np.char,('nchar'))
    methVar = ncid.createVariable('method','S1',('nchar'))

    aVar = ncid.createVariable('a',np.float32)
    bVar = ncid.createVariable('b',np.float32)
    ttVar = ncid.createVariable('tune_thresh',np.float32)
    srVar = ncid.createVariable('sm_rad',np.float32)
    fdVar = ncid.createVariable('fill_dbz',np.float32)
    bdVar = ncid.createVariable('bg_diff',np.float32)
    rtVar = ncid.createVariable('rain_type',np.int32,('time','lat','lon'),zlib=True,fill_value=missing_value )
    rtbVar = ncid.createVariable('rain_type_basic',np.int32,('time','lat','lon'),zlib=True,fill_value=missing_value )
    ccVar = ncid.createVariable('conv_cores',np.int32,('time','lat','lon'),
                                zlib=True,fill_value=missing_value )
    bgVar = ncid.createVariable('bkgrnd',np.float32,('time','lat','lon'),
                                zlib=True,fill_value=missing_value )

    # create variable attributes
    timeVar.standard_name = 'time'
    timeVar.long_name = 'Data time'
    timeVar.units = 'seconds since 1970-01-01T00:00:00Z'
    timeVar.calendar = 'standard'
    timeVar.comment = datetime
    
    lonVar.standard_name = 'longitude'
    lonVar.long_name = 'longitude'
    lonVar.units = 'degrees_east'
    
    latVar.standard_name = 'latitude'
    latVar.long_name = 'latitude'
    latVar.units = 'degrees_north'
        
    # ctVar
    ctVar.units = 'dBZ'
    ctVar.long_name = 'core_threshold'

    # methVar
    methVar.units = 'none'
    methVar.long_name = 'convstra_algorithm_used'
    methVar.comment = 'Either SYH (Steiner/Yuter/Houze) or YH (Yuter/Houze)'

    # aVar
    aVar.units = 'dBZ'
    aVar.long_name = 'tuning_variable1'
    aVar.comment = 'Tuning variable from YH [8. in YH]'

    # bVar
    bVar.units = 'dBZ'
    bVar.long_name = 'tuning variable2'
    bVar.comment = 'Tuning variable from YH [64. in YH]'

    # ttVar
    ttVar.units = 'dBZ'
    ttVar.long_name = 'tuning_threshold'
    ttVar.comment = 'tuning threshold'

    # srVar
    srVar.units = 'km'
    srVar.long_name = 'sm_radius'
    srVar.comment = 'sm_radius'

    # fdVar
    fdVar.units = 'dBZ'
    fdVar.long_name = 'fill_dbz'
    fdVar.comment = 'fill_dbz'

    # bdVar
    bdVar.units = 'dBZ'
    bdVar.long_name = 'shallow_conv_min'
    bdVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # rain_type
    rtVar.units = 'none'
    rtVar.long_name = 'rain_type_classification'
    rtVar.NO_ECHO = '0'
    rtVar.STRATIFORM = '1'
    rtVar.CONVECTIVE = '2'
    rtVar.UNCERTAIN = '3'
    rtVar.ancillary_variables = 'core_thresh method a b tune_thres sm_rad fill_dbz bg_diff'

    # rain_type basic
    rtbVar.units = 'none'
    rtbVar.long_name = 'rain_type_simple_classification'
    rtbVar.NO_ECHO = '0'
    rtbVar.STRATIFORM = '1'
    rtbVar.CONVECTIVE = '2'
    rtbVar.UNCERTAIN = '3'
    rtbVar.ancillary_variables = 'core_thresh method a b tune_thres sm_rad fill_dbz bg_diff'

    # conv_cores
    ccVar.units = 'none'
    ccVar.long_name = 'convective_cores'

    # bkgrnd
    bgVar.units = 'dBZ'
    bgVar.long_name = 'background_values'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.radar_latitude = lat_origin
    ncid.radar_longitude = lon_origin
    ncid.institution = institution
    ncid.source = source
    ncid.title = title
    ncid.references = references
    ncid.comment = comment
    ncid.history = 'File created ' + currentTime

    # write vars to file
    timeVar[:] = time_val
    lonVar[:] = lon_vals
    latVar[:] = lat_vals
    ctVar[:] = core_thresh
    methVar[:] = method
    aVar[:] = a
    bVar[:] = b
    ttVar[:] = tune_thresh
    srVar[:] = sm_rad
    fdVar[:] = fill_dbz
    bdVar[:] = bg_diff
    rtVar[:] = rain_type_vals[np.newaxis,:]
    rtbVar[:] = rain_type_basic_vals[np.newaxis,:]
    ccVar[:] = conv_core_vals[np.newaxis,:]
    bgVar[:] = bkgrnd_vals[np.newaxis,:]

    # close file
    ncid.close()
