import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

def writeBasic_RRNetcdf(ncname,aa_conv,bb_conv,aa_stra,bb_stra,
                        aa_blen,bb_blen,title,institution,source,
                        references,comment,rr_stra2uncert_name,
                        rr_conv2uncert_name,rr_blen2uncert_name,
                        rr_blen2uncert_stra2frin_name,dx,radar_lat,
                        radar_lon,xdim,ydim,rr_stra2uncert,
                        rr_conv2uncert,rr_blen2uncert,
                        rr_blen2uncert_stra2frin,missing_value):

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
    aaStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bbStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aaConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bbConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aaBlenVar = ncid.createVariable('rr_blen_coef',np.float32,zlib=True )
    bbBlenVar = ncid.createVariable('rr_blen_exp',np.float32,zlib=True )
    rrStra2uncertVar = ncid.createVariable(rr_stra2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrConv2uncertVar = ncid.createVariable(rr_conv2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBlen2uncertVar = ncid.createVariable(rr_blen2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBlen2uncertStra2frinVar = ncid.createVariable(rr_blen2uncert_stra2frin_name, \
                                                    np.float32,('time','y','x'),zlib=True, \
                                                    fill_value=missing_value )

    # create variable attributes
    xspVar.units = 'km'
    yspVar.units = 'km'

    aaStraVar.units = 'none'
    aaStraVar.long_name = 'coefficient in stratiform R-Z equation'
    bbStraVar.units = 'none'
    bbStraVar.long_name = 'exponent in stratiform R-Z equation'

    aaConvVar.units = 'none'
    aaConvVar.long_name = 'coefficient in convective R-Z equation'
    bbConvVar.units = 'none'
    bbConvVar.long_name = 'exponent in convective R-Z equation'

    aaBlenVar.units = 'none'
    aaBlenVar.long_name = 'coefficient in blended R-Z equation'
    bbBlenVar.units = 'none'
    bbBlenVar.long_name = 'exponent in blended R-Z equation'

    # rainrates
    rrStra2uncertVar.units = 'mm/hr'
    rrStra2uncertVar.long_name = 'rain_rate_stra2uncert'
    rrStra2uncertVar.coordinates = 'lon lat'
    rrStra2uncertVar.grid_mapping = 'grid_mapping'
    rrStra2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrConv2uncertVar.units = 'mm/hr'
    rrConv2uncertVar.long_name = 'rain_rate_conv2uncert'
    rrConv2uncertVar.coordinates = 'lon lat'
    rrConv2uncertVar.grid_mapping = 'grid_mapping'
    rrConv2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrBlen2uncertVar.units = 'mm/hr'
    rrBlen2uncertVar.long_name = 'rain_rate_blen2uncert'
    rrBlen2uncertVar.coordinates = 'lon lat'
    rrBlen2uncertVar.grid_mapping = 'grid_mapping'
    rrBlen2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_blen_coef rr_blen_exp'

    rrBlen2uncertStra2frinVar.units = 'mm/hr'
    rrBlen2uncertStra2frinVar.long_name = 'rain_rate_blen2uncert_stra2frin'
    rrBlen2uncertStra2frinVar.coordinates = 'lon lat'
    rrBlen2uncertStra2frinVar.grid_mapping = 'grid_mapping'
    rrBlen2uncertStra2frinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_blen_coef rr_blen_exp'

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
    aaStraVar[:] = aa_stra
    bbStraVar[:] = bb_stra
    aaConvVar[:] = aa_conv
    bbConvVar[:] = bb_conv
    aaBlenVar[:] = aa_blen
    bbBlenVar[:] = bb_blen
    rrStra2uncertVar[0,:,:] = rr_stra2uncert
    rrConv2uncertVar[0,:,:] = rr_conv2uncert
    rrBlen2uncertVar[0,:,:] = rr_blen2uncert
    rrBlen2uncertStra2frinVar[0,:,:] = rr_blen2uncert_stra2frin

    # close file
    ncid.close()

def writeCF_RRnetcdf(ncname,aa_conv,bb_conv,aa_stra,bb_stra,
                     aa_blen,bb_blen,title,institution,source,
                     references,comment,rr_stra2uncert_name,
                     rr_conv2uncert_name,rr_blen2uncert_name,
                     rr_blen2uncert_stra2frin_name,timeVal,xVal,
                     yVal,latVal,lonVal,gmVal,lat_origin,
                     lon_origin,rr_stra2uncert,rr_conv2uncert,
                     rr_blen2uncert,rr_blen2uncert_stra2frin,
                     missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_stra2uncert.shape[1])
    x = ncid.createDimension('x',rr_stra2uncert.shape[0])

    # create variables
    timeVar = ncid.createVariable('time',np.float64,('time'),zlib=True )
    xVar = ncid.createVariable('x',np.float32,('x'),zlib=True )
    yVar = ncid.createVariable('y',np.float32,('y'),zlib=True )
    latVar = ncid.createVariable('lat',np.float32,('y','x'),zlib=True )
    lonVar = ncid.createVariable('lon',np.float32,('y','x'),zlib=True )
    gmVar = ncid.createVariable('grid_mapping',np.int32)
    aaStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bbStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aaConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bbConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aaBlenVar = ncid.createVariable('rr_blen_coef',np.float32,zlib=True )
    bbBlenVar = ncid.createVariable('rr_blen_exp',np.float32,zlib=True )
    rrStra2uncertVar = ncid.createVariable(rr_stra2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrConv2uncertVar = ncid.createVariable(rr_conv2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBlen2uncertVar = ncid.createVariable(rr_blen2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBlen2uncertStra2frinVar = ncid.createVariable(rr_blen2uncert_stra2frin_name,np.float32, \
                                                    ('time','y','x'),zlib=True, \
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
    
    aaStraVar.units = 'none'
    aaStraVar.long_name = 'coefficient in stratiform R-Z equation'
    bbStraVar.units = 'none'
    bbStraVar.long_name = 'exponent in stratiform R-Z equation'

    aaConvVar.units = 'none'
    aaConvVar.long_name = 'coefficient in convective R-Z equation'
    bbConvVar.units = 'none'
    bbConvVar.long_name = 'exponent in convective R-Z equation'

    aaBlenVar.units = 'none'
    aaBlenVar.long_name = 'coefficient in blended R-Z equation'
    bbBlenVar.units = 'none'
    bbBlenVar.long_name = 'exponent in blended R-Z equation'

    # rainrates
    rrStra2uncertVar.units = 'mm/hr'
    rrStra2uncertVar.long_name = 'rain_rate_stra2uncert'
    rrStra2uncertVar.coordinates = 'lon lat'
    rrStra2uncertVar.grid_mapping = 'grid_mapping'
    rrStra2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrConv2uncertVar.units = 'mm/hr'
    rrConv2uncertVar.long_name = 'rain_rate_conv2uncert'
    rrConv2uncertVar.coordinates = 'lon lat'
    rrConv2uncertVar.grid_mapping = 'grid_mapping'
    rrConv2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrBlen2uncertVar.units = 'mm/hr'
    rrBlen2uncertVar.long_name = 'rain_rate_blen2uncert'
    rrBlen2uncertVar.coordinates = 'lon lat'
    rrBlen2uncertVar.grid_mapping = 'grid_mapping'
    rrBlen2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_blen_coef rr_blen_exp'

    rrBlen2uncertStra2frinVar.units = 'mm/hr'
    rrBlen2uncertStra2frinVar.long_name = 'rain_rate_blen2uncert_stra2frin'
    rrBlen2uncertStra2frinVar.coordinates = 'lon lat'
    rrBlen2uncertStra2frinVar.grid_mapping = 'grid_mapping'
    rrBlen2uncertStra2frinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_blen_coef rr_blen_exp'

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
    aaStraVar[:] = aa_stra
    bbStraVar[:] = bb_stra
    aaConvVar[:] = aa_conv
    bbConvVar[:] = bb_conv
    aaBlenVar[:] = aa_blen
    bbBlenVar[:] = bb_blen
    rrStra2uncertVar[0,:,:] = rr_stra2uncert
    rrConv2uncertVar[0,:,:] = rr_conv2uncert
    rrBlen2uncertVar[0,:,:] = rr_blen2uncert
    rrBlen2uncertStra2frinVar[0,:,:] = rr_blen2uncert_stra2frin

    # close file
    ncid.close()

def writeZeb_RRNetcdf(ncname,aa_conv,bb_conv,aa_stra,bb_stra,
                      aa_blen,bb_blen,title,institution,source,
                      references,comment,rr_stra2uncert_name,
                      rr_conv2uncert_name,rr_blen2uncert_name,
                      rr_blen2uncert_stra2frin_name,btVal,toVal,
                      latVal,lonVal,altVal,xspVal,yspVal,zspVal,
                      rr_stra2uncert,rr_conv2uncert,
                      rr_blen2uncert,rr_blen2uncert_stra2frin,
                      missing_value):

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
    aaStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bbStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aaConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bbConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aaBlenVar = ncid.createVariable('rr_blen_coef',np.float32,zlib=True )
    bbBlenVar = ncid.createVariable('rr_blen_exp',np.float32,zlib=True )
    rrStra2uncertVar = ncid.createVariable(rr_stra2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrConv2uncertVar = ncid.createVariable(rr_conv2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBlen2uncertVar = ncid.createVariable(rr_blen2uncert_name,np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBlen2uncertStra2frinVar = ncid.createVariable(rr_blen2uncert_stra2frin_name, \
                                                    np.float32,('time','y','x'),zlib=True, \
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
    
    aaStraVar.units = 'none'
    aaStraVar.long_name = 'coefficient in stratiform R-Z equation'
    bbStraVar.units = 'none'
    bbStraVar.long_name = 'exponent in stratiform R-Z equation'

    aaConvVar.units = 'none'
    aaConvVar.long_name = 'coefficient in convective R-Z equation'
    bbConvVar.units = 'none'
    bbConvVar.long_name = 'exponent in convective R-Z equation'

    aaBlenVar.units = 'none'
    aaBlenVar.long_name = 'coefficient in blended R-Z equation'
    bbBlenVar.units = 'none'
    bbBlenVar.long_name = 'exponent in blended R-Z equation'

    # rainrates
    rrStra2uncertVar.units = 'mm/hr'
    rrStra2uncertVar.long_name = 'rain_rate_stra2uncert'
    rrStra2uncertVar.coordinates = 'lon lat'
    rrStra2uncertVar.grid_mapping = 'grid_mapping'
    rrStra2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrConv2uncertVar.units = 'mm/hr'
    rrConv2uncertVar.long_name = 'rain_rate_conv2uncert'
    rrConv2uncertVar.coordinates = 'lon lat'
    rrConv2uncertVar.grid_mapping = 'grid_mapping'
    rrConv2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrBlen2uncertVar.units = 'mm/hr'
    rrBlen2uncertVar.long_name = 'rain_rate_blen2uncert'
    rrBlen2uncertVar.coordinates = 'lon lat'
    rrBlen2uncertVar.grid_mapping = 'grid_mapping'
    rrBlen2uncertVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_blen_coef rr_blen_exp'

    rrBlen2uncertStra2frinVar.units = 'mm/hr'
    rrBlen2uncertStra2frinVar.long_name = 'rain_rate_blen2uncert_stra2frin'
    rrBlen2uncertStra2frinVar.coordinates = 'lon lat'
    rrBlen2uncertStra2frinVar.grid_mapping = 'grid_mapping'
    rrBlen2uncertStra2frinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_blen_coef rr_blen_exp'

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
    aaStraVar[:] = aa_stra
    bbStraVar[:] = bb_stra
    aaConvVar[:] = aa_conv
    bbConvVar[:] = bb_conv
    aaBlenVar[:] = aa_blen
    bbBlenVar[:] = bb_blen
    rrStra2uncertVar[0,:,:] = rr_stra2uncert
    rrConv2uncertVar[0,:,:] = rr_conv2uncert
    rrBlen2uncertVar[0,:,:] = rr_blen2uncert
    rrBlen2uncertStra2frinVar[0,:,:] = rr_blen2uncert_stra2frin

    #close file
    ncid.close()
