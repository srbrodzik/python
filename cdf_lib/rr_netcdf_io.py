import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

def writeBasic_RRNetcdf(ncname,a_conv,b_conv,a_stra,b_stra,title,
                        institution,source,references,comment,dx,
                        radar_lat,radar_lon,xdim,ydim,rainrate_min,
                        rainrate_max,rr_2zr_min,rr_2zr_max,rr_1zr_min,
                        rr_1zr_max,rainrate_csu,method_csu,
                        missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rainrate_min.shape[1])
    x = ncid.createDimension('x',rainrate_min.shape[0])

    # create variables
    xspVar = ncid.createVariable('x_spacing',np.float32,zlib=True )
    yspVar = ncid.createVariable('y_spacing',np.float32,zlib=True )
    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )

    # create variable attributes

    xspVar.units = 'km'
    yspVar.units = 'km'
    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'rain_rate_minimum'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'
    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'rain_rate_maximum'
    rrMaxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    # create global attributes
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source = source
    ncid.references = references
    ncid.comment = comment
    ncid.radar_lat = radar_lat
    ncid.radar_lon = radar_lon

    # write variables to file
    xspVar[:] = dx
    yspVar[:] = dx
    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    rrMinVar[0,:,:] = rainrate_min
    rrMaxVar[0,:,:] = rainrate_max

    # close file
    ncid.close()

def writeCF_RRnetcdf(ncname,a_conv,b_conv,a_stra,b_stra,title,
                     institution,source,references,comment,
                     timeVal,xVal,yVal,latVal,lonVal,gmVal,
                     lat_origin,lon_origin,rainrate_min,
                     rainrate_max,rr_2zr_min,rr_2zr_max,
                     rr_1zr_min,rr_1zr_max,rainrate_csu,
                     method_csu,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rainrate_min.shape[1])
    x = ncid.createDimension('x',rainrate_min.shape[0])

    # create variables
    timeVar = ncid.createVariable('time',np.float64,('time'),zlib=True )
    xVar = ncid.createVariable('x',np.float32,('x'),zlib=True )
    yVar = ncid.createVariable('y',np.float32,('y'),zlib=True )
    latVar = ncid.createVariable('lat',np.float32,('y','x'),zlib=True )
    lonVar = ncid.createVariable('lon',np.float32,('y','x'),zlib=True )
    gmVar = ncid.createVariable('grid_mapping',np.int32)
    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rr1minVar = ncid.createVariable('rain_rate_1zr_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rr1maxVar = ncid.createVariable('rain_rate_1zr_max',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rr2minVar = ncid.createVariable('rain_rate_2zr_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rr2maxVar = ncid.createVariable('rain_rate_2zr_max',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrCsuVar = ncid.createVariable('rain_rate_csu',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    meCsuVar = ncid.createVariable('method_csu',np.int32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )

    # create variable attributes
    timeVar.standard_name = 'time'
    timeVar.long_name = 'Data time'
    timeVar.units = 'seconds since 1970-01-01T00:00:00Z'
    timeVar.calendar = 'standard'
    timeVar.axis = 'T'
    #timeVar.bounds = 'time_bounds'
    timeVar.comment = datetime
    
    xVar.standard_name = 'projection_x_coordinate'
    xVar.long_name = 'x distance on the projection plane from the origin'
    xVar.units = 'km'
    xVar.axis = 'X'
    
    yVar.standard_name = 'projection_y_coordinate'
    yVar.long_name = 'y distance on the projection plane from the origin'
    yVar.units = 'km'
    yVar.axis = 'Y'
    
    latVar.standard_name = 'latitude'
    latVar.units = 'degrees_north'
    
    lonVar.standard_name = 'longitude'
    lonVar.units = 'degrees_east'
    
    gmVar.grid_mapping_name = 'azimuthal_equidistant'
    gmVar.longitude_of_projection_origin = lon_origin
    gmVar.latitude_of_projection_origin = lat_origin
    gmVar.false_easting = 0
    gmVar.false_northing = 0
    
    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'

    # rainrates
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'rain_rate_minimum'
    #rrMinVar.grid_mapping = 'azimuthal_equidistant'
    rrMinVar.coordinates = 'lon lat'
    rrMinVar.grid_mapping = 'grid_mapping'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'rain_rate_maximum'
    #rrMaxVar.grid_mapping = 'azimuthal_equidistant'
    rrMaxVar.coordinates = 'lon lat'
    rrMaxVar.grid_mapping = 'grid_mapping'
    rrMaxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rr1minVar.units = 'mm/hr'
    rr1minVar.long_name = 'rain_rate_1zr_min'
    #rr1minVar.grid_mapping = 'azimuthal_equidistant'
    rr1minVar.coordinates = 'lon lat'
    rr1minVar.grid_mapping = 'grid_mapping'
    rr1minVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rr1maxVar.units = 'mm/hr'
    rr1maxVar.long_name = 'rain_rate_1zr_max'
    #rr1maxVar.grid_mapping = 'azimuthal_equidistant'
    rr1maxVar.coordinates = 'lon lat'
    rr1maxVar.grid_mapping = 'grid_mapping'
    rr1maxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rr2minVar.units = 'mm/hr'
    rr2minVar.long_name = 'rain_rate_2zr_min'
    #rr2minVar.grid_mapping = 'azimuthal_equidistant'
    rr2minVar.coordinates = 'lon lat'
    rr2minVar.grid_mapping = 'grid_mapping'
    rr2minVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rr2maxVar.units = 'mm/hr'
    rr2maxVar.long_name = 'rain_rate_2zr_max'
    #rr2maxVar.grid_mapping = 'azimuthal_equidistant'
    rr2maxVar.coordinates = 'lon lat'
    rr2maxVar.grid_mapping = 'grid_mapping'
    rr2maxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrCsuVar.units = 'mm/hr'
    rrCsuVar.long_name = 'rain_rate_csu_blended'
    #rrCsuVar.grid_mapping = 'azimuthal_equidistant'
    rrCsuVar.coordinates = 'lon lat'
    rrCsuVar.grid_mapping = 'grid_mapping'
    #rrCsuVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    meCsuVar.units = 'none'
    meCsuVar.long_name = 'method_csu_blended'
    #meCsuVar.grid_mapping = 'azimuthal_equidistant'
    meCsuVar.coordinates = 'lon lat'
    meCsuVar.grid_mapping = 'grid_mapping'
    meCsuVar.comment = '1=kdp_zdr, 2=kdp, 3=dz_zdr, 4=nexrad, 5=dz_rainonly'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source = source
    ncid.references = references
    ncid.comment = comment

    # write vars to file
    timeVar[:] = timeVal
    xVar[:] = xVal
    yVar[:] = yVal
    latVar[:] = latVal
    lonVar[:] = lonVal
    #gmVar[:] = gmVal
    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    rrMinVar[0,:,:] = rainrate_min
    rrMaxVar[0,:,:] = rainrate_max
    rr1minVar[0,:,:] = rr_1zr_min
    rr1maxVar[0,:,:] = rr_1zr_max
    rr2minVar[0,:,:] = rr_2zr_min
    rr2maxVar[0,:,:] = rr_2zr_max
    rrCsuVar[0,:,:] = rainrate_csu
    meCsuVar[0,:,:] = method_csu

    # close file
    ncid.close()

def writeZeb_RRNetcdf(ncname,a_conv,b_conv,a_stra,b_stra,title,
                      institution,source,references,comment,btVal,
                      toVal,latVal,lonVal,altVal,xspVal,yspVal,
                      zspVal,rainrate_min,rainrate_max,rr_2zr_min,
                      rr_2zr_max,rr_1zr_min,rr_1zr_max,rainrate_csu,
                      method_csu,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    z = ncid.createDimension('z',1)
    y = ncid.createDimension('y',rainrate_min.shape[1])
    x = ncid.createDimension('x',rainrate_min.shape[0])

    # create variables
    bt = ncid.createVariable('base_time',np.float64 )
    to = ncid.createVariable('time_offset',np.float32,('time',) )
    lat = ncid.createVariable('lat',np.float32 )
    lon = ncid.createVariable('lon',np.float32 )
    alt = ncid.createVariable('alt',np.float32 )
    xsp = ncid.createVariable('x_spacing',np.float32 )
    ysp = ncid.createVariable('y_spacing',np.float32 )
    zsp = ncid.createVariable('z_spacing',np.float32 )
    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','z','y','x',), \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','z','y','x',), \
                                   fill_value=missing_value )

    # create variable attributes
    bt.units = 'seconds since 1970-01-01 00:00:00 +0000'
    to.units = 'seconds since base_time'
    lat.units = 'degrees_north'
    lon.units = 'degrees_east'
    alt.units = 'km'
    xsp.units = 'km'
    ysp.units = 'km'
    zsp.units = 'km'
    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'rain_rate_minimum'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'
    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'rain_rate_maximum'
    rrMaxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    # create global attributes
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source = source
    ncid.references = references
    ncid.comment = comment

    # write vars to file
    bt[:] = btVal
    to[:] = toVal
    lat[:] = latVal
    lon[:] = lonVal
    alt[:] = altVal
    xsp[:] = xspVal
    ysp[:] = yspVal
    zsp[:] = zspVal
    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    rrMinVar[0,0,:,:] = rainrate_min
    rrMaxVar[0,0,:,:] = rainrate_max

    #close file
    ncid.close()
