import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

def writeBasicNetcdf(ncname,rt_min,rt_best,rt_max,rr_min,rr_best,rr_max,
                     method_min,method_best,method_max,
                     deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                     weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                     a_conv,b_conv,a_stra,b_stra,
                     title,institution,source_uw,source_csu,references,comment,
                     dx,radar_lat,radar_lon,xdim,ydim,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_min.shape[1])
    x = ncid.createDimension('x',rr_min.shape[0])

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
    rrMinVar.long_name = 'minimum_uncertainity_value_of_radar_estimated_rainrate'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'
    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'maximum_uncertainity_value_of_radar_estimated_rainrate'
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
    rrMinVar[0,:,:] = rr_min
    rrMaxVar[0,:,:] = rr_max

    # close file
    ncid.close()

def writeBasicNetcdfV6(ncname,rt,rr_min,rr_best,rr_max,
                       method_min,method_best,method_max,
                       deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                       weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                       a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                       a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                       title,institution,source_uw,source_csu,references,comment,
                       dx,radar_lat,radar_lon,xdim,ydim,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_min.shape[1])
    x = ncid.createDimension('x',rr_min.shape[0])

    # create variables
    xspVar = ncid.createVariable('x_spacing',np.float32,zlib=True )
    yspVar = ncid.createVariable('y_spacing',np.float32,zlib=True )

    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )

    rtVar = ncid.createVariable('rain_type',np.int32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )

    # create variable attributes

    xspVar.units = 'km'
    yspVar.units = 'km'
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'

    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintypes
    rtVar.units = 'none'
    rtVar.long_name = 'rain_type_classification'
    rtVar.coordinates = 'lon lat'
    rtVar.grid_mapping = 'grid_mapping'
    rtVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                     types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                     types['WEAK_ECHO']))
    rtVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrates
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'minimum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMinVar.grid_mapping = 'azimuthal_equidistant'
    rrMinVar.coordinates = 'lon lat'
    rrMinVar.grid_mapping = 'grid_mapping'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    #rrBestVar.grid_mapping = 'azimuthal_equidistant'
    rrBestVar.coordinates = 'lon lat'
    rrBestVar.grid_mapping = 'grid_mapping'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'maximum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMaxVar.grid_mapping = 'azimuthal_equidistant'
    rrMaxVar.coordinates = 'lon lat'
    rrMaxVar.grid_mapping = 'grid_mapping'
    rrMaxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

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
    
    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rtVar[0,:,:] = rt

    rrMinVar[0,:,:] = rr_min
    rrBestVar[0,:,:] = rr_best
    rrMaxVar[0,:,:] = rr_max

    # close file
    ncid.close()

def writeCFnetcdf(ncname,rt_min,rt_best,rt_max,rr_min,rr_best,rr_max,
                  method_min,method_best,method_max,
                  deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                  weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                  a_conv,b_conv,a_stra,b_stra,
                  title,institution,source_uw,source_csu,references,comment,
                  timeVal,xVal,yVal,latVal,lonVal,gmVal,
                  lat_origin,lon_origin,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_min.shape[1])
    x = ncid.createDimension('x',rr_min.shape[0])

    # create variables
    timeVar = ncid.createVariable('time',np.float64,('time'),zlib=True )
    xVar = ncid.createVariable('x',np.float32,('x'),zlib=True )
    yVar = ncid.createVariable('y',np.float32,('y'),zlib=True )
    latVar = ncid.createVariable('lat',np.float32,('y','x'),zlib=True )
    lonVar = ncid.createVariable('lon',np.float32,('y','x'),zlib=True )
    gmVar = ncid.createVariable('grid_mapping',np.int32)
    
    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    
    rtMinVar = ncid.createVariable('rain_type_min',np.int32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rtBestVar = ncid.createVariable('rain_type',np.int32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rtMaxVar = ncid.createVariable('rain_type_max',np.int32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','y','x'),zlib=True, \
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
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'

    # raintypes
    rtMinVar.units = 'none'
    rtMinVar.long_name = 'rain_type_classification_minimum'
    #rtMinVar.grid_mapping = 'azimuthal_equidistant'
    rtMinVar.coordinates = 'lon lat'
    rtMinVar.grid_mapping = 'grid_mapping'
    rtMinVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                     types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                     types['WEAK_ECHO']))
    rtMinVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtMinVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    rtBestVar.units = 'none'
    rtBestVar.long_name = 'rain_type_classification'
    #rtBestVar.grid_mapping = 'azimuthal_equidistant'
    rtBestVar.coordinates = 'lon lat'
    rtBestVar.grid_mapping = 'grid_mapping'
    rtBestVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                     types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                     types['WEAK_ECHO']))
    rtBestVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtBestVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    rtMaxVar.units = 'none'
    rtMaxVar.long_name = 'rain_type_classification_maximum'
    #rtMaxVar.grid_mapping = 'azimuthal_equidistant'
    rtMaxVar.coordinates = 'lon lat'
    rtMaxVar.grid_mapping = 'grid_mapping'
    rtMaxVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                     types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                     types['WEAK_ECHO']))
    rtMaxVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtMaxVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'
    
    # rainrates
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'minimum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMinVar.grid_mapping = 'azimuthal_equidistant'
    rrMinVar.coordinates = 'lon lat'
    rrMinVar.grid_mapping = 'grid_mapping'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    #rrBestVar.grid_mapping = 'azimuthal_equidistant'
    rrBestVar.coordinates = 'lon lat'
    rrBestVar.grid_mapping = 'grid_mapping'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'maximum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMaxVar.grid_mapping = 'azimuthal_equidistant'
    rrMaxVar.coordinates = 'lon lat'
    rrMaxVar.grid_mapping = 'grid_mapping'
    rrMaxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
    ncid.source_csu = source_csu
    ncid.references = references
    ncid.comment = comment

    # write vars to file
    timeVar[:] = timeVal
    xVar[:] = xVal
    yVar[:] = yVal
    latVar[:] = latVal
    lonVar[:] = lonVal
    #gmVar[:] = gmVal

    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra

    rtMinVar[0,:,:] = rt_min
    rtBestVar[0,:,:] = rt_best
    rtMaxVar[0,:,:] = rt_max

    rrMinVar[0,:,:] = rr_min
    rrBestVar[0,:,:] = rr_best
    rrMaxVar[0,:,:] = rr_max

    # close file
    ncid.close()

def writeCFnetcdfV6(ncname,rt,rr_min,rr_best,rr_max,
                    method_min,method_best,method_max,zr_method,
                    deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                    weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                    a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                    a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                    title,institution,source_uw,source_csu,references,comment,
                    timeVal,xVal,yVal,latVal,lonVal,gmVal,
                    lat_origin,lon_origin,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_min.shape[1])
    x = ncid.createDimension('x',rr_min.shape[0])

    # create variables
    timeVar = ncid.createVariable('time',np.float64,('time'),zlib=True )
    xVar = ncid.createVariable('x',np.float32,('x'),zlib=True )
    yVar = ncid.createVariable('y',np.float32,('y'),zlib=True )
    latVar = ncid.createVariable('lat',np.float32,('y','x'),zlib=True )
    lonVar = ncid.createVariable('lon',np.float32,('y','x'),zlib=True )
    gmVar = ncid.createVariable('grid_mapping',np.int32)
    
    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )
   
    rtVar = ncid.createVariable('rain_type',np.int32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    methMinVar = ncid.createVariable('zr_method_min',np.float32,('time','y','x'),zlib=True, \
                                     fill_value=-1 )
    methBestVar = ncid.createVariable('zr_method',np.float32,('time','y','x'),zlib=True, \
                                      fill_value=-1 )
    methMaxVar = ncid.createVariable('zr_method_max',np.float32,('time','y','x'),zlib=True, \
                                     fill_value=-1 )

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
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'

    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'

    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintype
    rtVar.units = 'none'
    rtVar.long_name = 'rain_type_classification'
    #rtVar.grid_mapping = 'azimuthal_equidistant'
    rtVar.coordinates = 'lon lat'
    rtVar.grid_mapping = 'grid_mapping'
    rtVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                     types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                     types['WEAK_ECHO']))
    rtVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrates
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'minimum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMinVar.grid_mapping = 'azimuthal_equidistant'
    rrMinVar.coordinates = 'lon lat'
    rrMinVar.grid_mapping = 'grid_mapping'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    #rrBestVar.grid_mapping = 'azimuthal_equidistant'
    rrBestVar.coordinates = 'lon lat'
    rrBestVar.grid_mapping = 'grid_mapping'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'maximum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMaxVar.grid_mapping = 'azimuthal_equidistant'
    rrMaxVar.coordinates = 'lon lat'
    rrMaxVar.grid_mapping = 'grid_mapping'
    rrMaxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    # zr methods
    methMinVar.units = ''
    methMinVar.long_name = 'zr_method_used_for_minimum_rainrate'
    #methMinVar.grid_mapping = 'azimuthal_equidistant'
    methMinVar.coordinates = 'lon lat'
    methMinVar.grid_mapping = 'grid_mapping'
    methMinVar.flag_values = np.array((zr_method['kdp_zdr'],zr_method['kdp'],zr_method['z_zdr'],
                                       zr_method['z_all'],zr_method['z_conv'],zr_method['z_stra']))
    methMinVar.flag_meanings = "kdp_zdr kdp z_zdr z_all z_conv z_stra"

    methBestVar.units = ''
    methBestVar.long_name = 'zr_method_used_for_rainrate'
    #methBestVar.grid_mapping = 'azimuthal_equidistant'
    methBestVar.coordinates = 'lon lat'
    methBestVar.grid_mapping = 'grid_mapping'
    methBestVar.flag_values = np.array((zr_method['kdp_zdr'],zr_method['kdp'],zr_method['z_zdr'],
                                        zr_method['z_all'],zr_method['z_conv'],zr_method['z_stra']))
    #methBestVar.flag_meanings = np.array(['kdp_zdr   kdp   z_zdr   z_all   z_conv   z_stra'])
    methBestVar.flag_meanings = "kdp_zdr kdp z_zdr z_all z_conv z_stra"

    methMaxVar.units = ''
    methMaxVar.long_name = 'zr_method_used_for_maximum_rainrate'
    #methMaxVar.grid_mapping = 'azimuthal_equidistant'
    methMaxVar.coordinates = 'lon lat'
    methMaxVar.grid_mapping = 'grid_mapping'
    methMaxVar.flag_values = np.array((zr_method['kdp_zdr'],zr_method['kdp'],zr_method['z_zdr'],
                                       zr_method['z_all'],zr_method['z_conv'],zr_method['z_stra']))
    methMaxVar.flag_meanings = "kdp_zdr kdp z_zdr z_all z_conv z_stra"

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
    ncid.source_csu = source_csu
    ncid.references = references
    ncid.comment = comment

    # write vars to file
    timeVar[:] = timeVal
    xVar[:] = xVal
    yVar[:] = yVal
    latVar[:] = latVal
    lonVar[:] = lonVal
    #gmVar[:] = gmVal

    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rtVar[0,:,:] = rt

    rrMinVar[0,:,:] = rr_min
    rrBestVar[0,:,:] = rr_best
    rrMaxVar[0,:,:] = rr_max

    methMinVar[0,:,:] = method_min
    methBestVar[0,:,:] = method_best
    methMaxVar[0,:,:] = method_max

    # close file
    ncid.close()

def writeCFnetcdfV6_rr_only(ncname,rr_min,rr_best,rr_max,
                            a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                            a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                            title,institution,source,references,comment,
                            timeVal,xVal,yVal,latVal,lonVal,gmVal,
                            lat_origin,lon_origin,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_min.shape[1])
    x = ncid.createDimension('x',rr_min.shape[0])

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
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )
   
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','y','x'),zlib=True, \
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

    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # rainrates
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'minimum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMinVar.grid_mapping = 'azimuthal_equidistant'
    rrMinVar.coordinates = 'lon lat'
    rrMinVar.grid_mapping = 'grid_mapping'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    #rrBestVar.grid_mapping = 'azimuthal_equidistant'
    rrBestVar.coordinates = 'lon lat'
    rrBestVar.grid_mapping = 'grid_mapping'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'maximum_uncertainity_value_of_radar_estimated_rainrate'
    #rrMaxVar.grid_mapping = 'azimuthal_equidistant'
    rrMaxVar.coordinates = 'lon lat'
    rrMaxVar.grid_mapping = 'grid_mapping'
    rrMaxVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.institution = institution
    ncid.source = source
    ncid.title = title
    ncid.references = references
    ncid.comment = comment
    ncid.history = 'File created ' + currentTime

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
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rrMinVar[0,:,:] = rr_min
    rrBestVar[0,:,:] = rr_best
    rrMaxVar[0,:,:] = rr_max

    # close file
    ncid.close()

def writeZebNetcdf(ncname,rt_min,rt_best,rt_max,rr_min,rr_best,rr_max,
                   method_min,method_best,method_max,
                   deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                   weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                   a_conv,b_conv,a_stra,b_stra,
                   title,institution,source_uw,source_csu,references,comment,
                   btVal,toVal,latVal,lonVal,altVal,xspVal,yspVal,zspVal,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    z = ncid.createDimension('z',1)
    y = ncid.createDimension('y',rr_min.shape[1])
    x = ncid.createDimension('x',rr_min.shape[0])

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
    rrMinVar.long_name = 'minimum_uncertainity_value_of_radar_estimated_rainrate'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'
    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'maximum_uncertainity_value_of_radar_estimated_rainrate'
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
    rrMinVar[0,0,:,:] = rr_min
    rrMaxVar[0,0,:,:] = rr_max

    #close file
    ncid.close()


def writeZebNetcdfV6(ncname,rt,rr_min,rr_best,rr_max,
                     method_min,method_best,method_max,
                     deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                     weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                     a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                     a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                     title,institution,source_uw,source_csu,references,comment,
                     btVal,toVal,latVal,lonVal,altVal,xspVal,yspVal,zspVal,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    z = ncid.createDimension('z',1)
    y = ncid.createDimension('y',rr_min.shape[1])
    x = ncid.createDimension('x',rr_min.shape[0])

    # create variables
    bt = ncid.createVariable('base_time',np.float64 )
    to = ncid.createVariable('time_offset',np.float32,('time',) )
    lat = ncid.createVariable('lat',np.float32 )
    lon = ncid.createVariable('lon',np.float32 )
    alt = ncid.createVariable('alt',np.float32 )
    xsp = ncid.createVariable('x_spacing',np.float32 )
    ysp = ncid.createVariable('y_spacing',np.float32 )
    zsp = ncid.createVariable('z_spacing',np.float32 )
    
    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )

    rtVar = ncid.createVariable('rain_type',np.int32,('time','z','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrMinVar = ncid.createVariable('rain_rate_min',np.float32,('time','z','y','x'), \
                                   fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','z','y','x'), \
                                   fill_value=missing_value )
    rrMaxVar = ncid.createVariable('rain_rate_max',np.float32,('time','z','y','x'), \
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
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'

    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintypes
    rtVar.units = 'none'
    rtVar.long_name = 'rain_type_classification'
    #rtVar.grid_mapping = 'azimuthal_equidistant'
    rtVar.coordinates = 'lon lat'
    rtVar.grid_mapping = 'grid_mapping'
    rtVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                     types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                     types['WEAK_ECHO']))
    rtVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrates
    rrMinVar.units = 'mm/hr'
    rrMinVar.long_name = 'minimum_uncertainity_value_of_radar_estimated_rainrate'
    rrMinVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp'

    rrMaxVar.units = 'mm/hr'
    rrMaxVar.long_name = 'maximum_uncertainity_value_of_radar_estimated_rainrate'
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
    
    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rtVar[0,0,:,:] = rt

    rrMinVar[0,0,:,:] = rr_min
    rrBestVar[0,0,:,:] = rr_best
    rrMaxVar[0,0,:,:] = rr_max

    #close file
    ncid.close()


def writeBasicNetcdfBestOnlyV6(ncname,rt,rr_best,method_best,
                               deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                               weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                               a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                               a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                               title,institution,source_uw,source_csu,references,comment,
                               dx,radar_lat,radar_lon,xdim,ydim,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_best.shape[1])
    x = ncid.createDimension('x',rr_best.shape[0])

    # create variables
    xspVar = ncid.createVariable('x_spacing',np.float32,zlib=True )
    yspVar = ncid.createVariable('y_spacing',np.float32,zlib=True )
    
    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )
    
    rtVar = ncid.createVariable('rain_type',np.float32,('time','y','x',),zlib=True, \
                                fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','y','x'),zlib=True, \
                                    fill_value=missing_value )

    # create variable attributes

    xspVar.units = 'km'
    yspVar.units = 'km'
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'
    
    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintype
    rtVar.units = 'none'
    rtVar.long_name = 'rain_type_classification'
    #rtVar.grid_mapping = 'azimuthal_equidistant'
    rtVar.coordinates = 'lon lat'
    rtVar.grid_mapping = 'grid_mapping'
    rtVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                  types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                  types['WEAK_ECHO']))
    rtVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrate
    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

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
    
    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp
    
    rtVar[0,:,:] = rt

    rrBestVar[0,:,:] = rr_best

    # close file
    ncid.close()


def writeCFnetcdfBestOnly(ncname,rt_best,rr_best,method_best,
                          deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                          weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                          a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                          a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                          title,institution,source_uw,source_csu,references,comment,
                          timeVal,xVal,yVal,latVal,lonVal,gmVal,
                          lat_origin,lon_origin,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_best.shape[1])
    x = ncid.createDimension('x',rr_best.shape[0])

    # create variables
    timeVar = ncid.createVariable('time',np.float64,('time'),zlib=True )
    xVar = ncid.createVariable('x',np.float32,('x'),zlib=True )
    yVar = ncid.createVariable('y',np.float32,('y'),zlib=True )
    latVar = ncid.createVariable('lat',np.float32,('y','x'),zlib=True )
    lonVar = ncid.createVariable('lon',np.float32,('y','x'),zlib=True )
    gmVar = ncid.createVariable('grid_mapping',np.int32)
    
    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )
    
    rtBestVar = ncid.createVariable('rain_type',np.int32,('time','y','x'),zlib=True, \
                                   fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','y','x'),zlib=True, \
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
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'

    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintype
    rtBestVar.units = 'none'
    rtBestVar.long_name = 'rain_type_classification'
    #rtBestVar.grid_mapping = 'azimuthal_equidistant'
    rtBestVar.coordinates = 'lon lat'
    rtBestVar.grid_mapping = 'grid_mapping'
    rtBestVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                     types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                     types['WEAK_ECHO']))
    rtBestVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtBestVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrate
    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    #rrBestVar.grid_mapping = 'azimuthal_equidistant'
    rrBestVar.coordinates = 'lon lat'
    rrBestVar.grid_mapping = 'grid_mapping'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
    ncid.source_csu = source_csu
    ncid.references = references
    ncid.comment = comment

    # write vars to file
    timeVar[:] = timeVal
    xVar[:] = xVal
    yVar[:] = yVal
    latVar[:] = latVal
    lonVar[:] = lonVal
    #gmVar[:] = gmVal

    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rtBestVar[0,:,:] = rt_best

    rrBestVar[0,:,:] = rr_best

    # close file
    ncid.close()

def writeCFnetcdfBestOnlyV6(ncname,rt_best,rr_best,method_best,
                            deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                            weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                            a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                            a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                            title,institution,source_uw,source_csu,references,comment,
                            timeVal,xVal,yVal,latVal,lonVal,gmVal,
                            lat_origin,lon_origin,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',rr_best.shape[1])
    x = ncid.createDimension('x',rr_best.shape[0])

    # create variables
    timeVar = ncid.createVariable('time',np.float64,('time'),zlib=True )
    xVar = ncid.createVariable('x',np.float32,('x'),zlib=True )
    yVar = ncid.createVariable('y',np.float32,('y'),zlib=True )
    latVar = ncid.createVariable('lat',np.float32,('y','x'),zlib=True )
    lonVar = ncid.createVariable('lon',np.float32,('y','x'),zlib=True )
    gmVar = ncid.createVariable('grid_mapping',np.int32)
    
    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )
    
    rtVar = ncid.createVariable('rain_type',np.int32,('time','y','x'),zlib=True, \
                                fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate',np.float32,('time','y','x'),zlib=True, \
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
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'
    
    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'

    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintype
    rtVar.units = 'none'
    rtVar.long_name = 'rain_type_classification'
    #rtVar.grid_mapping = 'azimuthal_equidistant'
    rtVar.coordinates = 'lon lat'
    rtVar.grid_mapping = 'grid_mapping'
    rtVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                  types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                  types['WEAK_ECHO']))
    rtVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrate
    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate'
    #rrBestVar.grid_mapping = 'azimuthal_equidistant'
    rrBestVar.coordinates = 'lon lat'
    rrBestVar.grid_mapping = 'grid_mapping'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
    ncid.source_csu = source_csu
    ncid.references = references
    ncid.comment = comment

    # write vars to file
    timeVar[:] = timeVal
    xVar[:] = xVal
    yVar[:] = yVal
    latVar[:] = latVal
    lonVar[:] = lonVal
    #gmVar[:] = gmVal

    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rtVar[0,:,:] = rt_best

    rrBestVar[0,:,:] = rr_best

    # close file
    ncid.close()

def writeZebNetcdfBestOnly(ncname,rt_best,rr_best,method_best,
                           deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                           weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                           a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                           a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                           title,institution,source_uw,source_csu,references,comment,
                           btVal,toVal,latVal,lonVal,altVal,xspVal,yspVal,zspVal,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    z = ncid.createDimension('z',1)
    y = ncid.createDimension('y',rr_best.shape[1])
    x = ncid.createDimension('x',rr_best.shape[0])

    # create variables
    bt = ncid.createVariable('base_time',np.float64 )
    to = ncid.createVariable('time_offset',np.float32,('time',) )
    lat = ncid.createVariable('lat',np.float32 )
    lon = ncid.createVariable('lon',np.float32 )
    alt = ncid.createVariable('alt',np.float32 )
    xsp = ncid.createVariable('x_spacing',np.float32 )
    ysp = ncid.createVariable('y_spacing',np.float32 )
    zsp = ncid.createVariable('z_spacing',np.float32 )

    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )

    rtBestVar = ncid.createVariable('rain_type_best',np.float32,('time','z','y','x',), \
                                    fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate_best',np.float32,('time','z','y','x',), \
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
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'

    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'
    
    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintype
    rtBestVar.units = 'none'
    rtBestVar.long_name = 'rain_type_classification_best'
    rtBestVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                      types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                      types['WEAK_ECHO']))
    rtBestVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtBestVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrate
    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate_best'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    # create global attributes
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
    ncid.source_csu = source_csu
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
    
    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rtBestVar[0,0,:,:] = rt_best

    rrBestVar[0,0,:,:] = rr_best

    #close file
    ncid.close()

def writeZebNetcdfBestOnlyV6(ncname,rt_best,rr_best,method_best,
                             deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                             weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                             a_conv,b_conv,a_stra,b_stra,a_all,b_all,a_kdp_zdr,b_KDP_zdr,b_kdp_ZDR,
                             a_zh_zdr,b_ZH_zdr,b_zh_ZDR,a_kdp,b_kdp,
                             title,institution,source_uw,source_csu,references,comment,
                             btVal,toVal,latVal,lonVal,altVal,xspVal,yspVal,zspVal,types,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    z = ncid.createDimension('z',1)
    y = ncid.createDimension('y',rr_best.shape[1])
    x = ncid.createDimension('x',rr_best.shape[0])

    # create variables
    bt = ncid.createVariable('base_time',np.float64 )
    to = ncid.createVariable('time_offset',np.float32,('time',) )
    lat = ncid.createVariable('lat',np.float32 )
    lon = ncid.createVariable('lon',np.float32 )
    alt = ncid.createVariable('alt',np.float32 )
    xsp = ncid.createVariable('x_spacing',np.float32 )
    ysp = ncid.createVariable('y_spacing',np.float32 )
    zsp = ncid.createVariable('z_spacing',np.float32 )

    zthVar = ncid.createVariable('rt_Z_th',np.float32)
    rbgVar = ncid.createVariable('rt_R_bg',np.float32)
    aVar = ncid.createVariable('rt_a',np.float32)
    bVar = ncid.createVariable('rt_b',np.float32)
    rcVar = ncid.createVariable('rt_R_conv',np.float32)
    zcVar = ncid.createVariable('rt_Z_conv',np.float32)
    zwVar = ncid.createVariable('rt_Z_weak',np.float32)
    zsVar = ncid.createVariable('rt_Z_shallow',np.float32)
    alVar = ncid.createVariable('rt_A_low',np.float32)
    amVar = ncid.createVariable('rt_A_med',np.float32)
    ahVar = ncid.createVariable('rt_A_high',np.float32)

    aConvVar = ncid.createVariable('rr_conv_coef',np.float32,zlib=True )
    bConvVar = ncid.createVariable('rr_conv_exp',np.float32,zlib=True )
    aStraVar = ncid.createVariable('rr_stra_coef',np.float32,zlib=True )
    bStraVar = ncid.createVariable('rr_stra_exp',np.float32,zlib=True )
    aAllVar = ncid.createVariable('rr_all_coef',np.float32,zlib=True )
    bAllVar = ncid.createVariable('rr_all_exp',np.float32,zlib=True )
    aKdpZdrVar = ncid.createVariable('rr_kdp_zdr_coef',np.float32,zlib=True )
    bKDPZdrVar = ncid.createVariable('rr_KDP_zdr_exp',np.float32,zlib=True )
    bKdpZDRVar = ncid.createVariable('rr_kdp_ZDR_exp',np.float32,zlib=True )
    aZhZdrVar = ncid.createVariable('rr_zh_zdr_coef',np.float32,zlib=True )
    bZHZdrVar = ncid.createVariable('rr_ZH_zdr_exp',np.float32,zlib=True )
    bZhZDRVar = ncid.createVariable('rr_zh_ZDR_exp',np.float32,zlib=True )
    aKdpVar = ncid.createVariable('rr_kdp_coef',np.float32,zlib=True )
    bKdpVar = ncid.createVariable('rr_kdp_exp',np.float32,zlib=True )

    rtBestVar = ncid.createVariable('rain_type_best',np.float32,('time','z','y','x',), \
                                    fill_value=missing_value )
    rrBestVar = ncid.createVariable('rain_rate_best',np.float32,('time','z','y','x',), \
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
    
    # Z_th
    zthVar.units = 'dBZ'
    zthVar.long_name = 'trunc_Z_conv_thres'
    zthVar.comment = 'reflectivity threshold at or above which echos are classified as convective'

    # R_bg
    rbgVar.units = 'km'
    rbgVar.long_name = 'backgrnd_radius'
    rbgVar.comment = 'radius within which background reflectivity is computed'

    # a
    aVar.units = 'dBZ'
    aVar.long_name = 'min_Z_diff'
    aVar.comment = 'factor for comparing echo to background reflectivity; see equation (1) in journal article referenced in "references" general attribute'

    # b
    bVar.units = 'dBZ'
    bVar.long_name = 'deep_cos_zero'
    bVar.comment = 'see equation (1) in journal article referenced in  "references" general attribute'

    # R_conv
    rcVar.units = 'km'
    rcVar.long_name = 'max_conv_radius'
    rcVar.comment = 'maximum radius around convective core for possible uncertain classification'

    # Z_conv
    zcVar.units = 'dBZ'
    zcVar.long_name = 'dbz_for_max_conv_radius'
    zcVar.comment = 'minimum dBZ required for max_conv_radius to apply'

    # Z_weak
    zwVar.units = 'dBZ'
    zwVar.long_name = 'weak_echo_thres'
    zwVar.comment = 'minimum dBZ for classification as not weak echo'

    # Z_shallow
    zsVar.units = 'dBZ'
    zsVar.long_name = 'shallow_conv_min'
    zsVar.comment = 'minimum dBZ for classification as convective for objects with area less than A-med'

    # A_low
    alVar.units = 'km^2'
    alVar.long_name = 'min_size'
    #alVar.comment = 'minimum areal coverage of echo object for classification as convective or stratiform'
    alVar.comment = 'minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV classification'

    # A_med
    amVar.units = 'km^2'
    amVar.long_name = 'start_slope'
    #amVar.comment = 'maximum areal coverage of echo object for allowing new Z_th equal to Z_shallow'
    amVar.comment = 'any contiguous echo object with areal coverage greater than this but less than A_high gets a new Z_th that is linearly interpolated between Z_shallow and Z_th depending on where area is between A_med and A_high'

    # A_high
    ahVar.units = 'km^2'
    ahVar.long_name = 'max_size'
    ahVar.comment = 'any contiguous echo object greater than this size gets a convective threshold of truncZconvthres'

    aConvVar.units = 'none'
    aConvVar.long_name = 'coefficient in convective Z-R equation'
    bConvVar.units = 'none'
    bConvVar.long_name = 'exponent in convective Z-R equation'

    aStraVar.units = 'none'
    aStraVar.long_name = 'coefficient in stratiform Z-R equation'
    bStraVar.units = 'none'
    bStraVar.long_name = 'exponent in stratiform Z-R equation'
    
    aAllVar.units = 'none'
    aAllVar.long_name = 'coefficient in single Z-R equation'
    bAllVar.units = 'none'
    bAllVar.long_name = 'exponent in single Z-R equation'

    aKdpZdrVar.units = 'none'
    aKdpZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bKDPZdrVar.units = 'none'
    bKDPZdrVar.long_name = 'kdp exponent in kdp-zdr Z-R equation'
    bKdpZDRVar.units = 'none'
    bKdpZDRVar.long_name = 'zdr exponent in kdp-zdr Z-R equation'

    aZhZdrVar.units = 'none'
    aZhZdrVar.long_name = 'coefficient in kdp-zdr Z-R equation'
    bZHZdrVar.units = 'none'
    bZHZdrVar.long_name = 'zh exponent in zh-zdr Z-R equation'
    bZhZDRVar.units = 'none'
    bZhZDRVar.long_name = 'zdr exponent in zh-zdr Z-R equation'

    aKdpVar.units = 'none'
    aKdpVar.long_name = 'coefficient in kdp Z-R equation'
    bKdpVar.units = 'none'
    bKdpVar.long_name = 'exponent in kdp Z-R equation'

    # raintype
    rtBestVar.units = 'none'
    rtBestVar.long_name = 'rain_type_classification_best'
    rtBestVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                      types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                      types['WEAK_ECHO']))
    rtBestVar.flag_meanings = "NO_SFC_ECHO STRATIFORM CONVECTIVE UNCERTAIN ISO_CONV_CORE ISO_CONV_FRINGE WEAK_ECHO"
    rtBestVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # rainrate
    rrBestVar.units = 'mm/hr'
    rrBestVar.long_name = 'radar_estimated_rainrate_best'
    rrBestVar.ancillary_variables = 'rr_conv_coef rr_conv_exp rr_stra_coef rr_stra_exp rr_all_coef rr_all_exp rr_kdp_zdr_coef rr_KDP_zdr_exp rr_kdp_ZDR_exp rr_zh_zdr_coef rr_ZH_zdr_exp rr_zh_ZDR_exp rr_kdp_coef rr_kdp_exp'

    # create global attributes
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
    ncid.source_csu = source_csu
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
    
    zthVar[:] = truncZconvthres
    rbgVar[:] = backgrndradius
    aVar[:] = minZdiff
    bVar[:] = deepcoszero
    rcVar[:] = maxConvRadius
    zcVar[:] = dBZformaxconvradius
    zwVar[:] = weakechothres
    zsVar[:] = shallowconvmin
    alVar[:] = minsize
    amVar[:] = startslope
    ahVar[:] = maxsize

    aConvVar[:] = a_conv
    bConvVar[:] = b_conv
    aStraVar[:] = a_stra
    bStraVar[:] = b_stra
    aAllVar[:] = a_all
    bAllVar[:] = b_all
    aKdpZdrVar[:] = a_kdp_zdr
    bKDPZdrVar[:] = b_KDP_zdr
    bKdpZDRVar[:] = b_kdp_ZDR
    aZhZdrVar[:] = a_zh_zdr
    bZHZdrVar[:] = b_ZH_zdr
    bZhZDRVar[:] = b_zh_ZDR
    aKdpVar[:] = a_kdp
    bKdpVar[:] = b_kdp

    rtBestVar[0,0,:,:] = rt_best

    rrBestVar[0,0,:,:] = rr_best

    #close file
    ncid.close()


