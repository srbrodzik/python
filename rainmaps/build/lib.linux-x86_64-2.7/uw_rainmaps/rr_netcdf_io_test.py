import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

def write_rtCFnetcdf(ncname,rt,
                     deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                     weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                     title,institution,source_uw,references,comment,
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
    y = ncid.createDimension('y',rt.shape[1])
    x = ncid.createDimension('x',rt.shape[0])

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

    rtVar = ncid.createVariable('rain_type',np.int32,('time','y','x'),zlib=True, \
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

    # raintype
    rtVar.units = 'none'
    rtVar.long_name = 'rain_type_classification'
    #rtVar.grid_mapping = 'azimuthal_equidistant'
    rtVar.coordinates = 'lon lat'
    rtVar.grid_mapping = 'grid_mapping'
    rtVar.flag_values = np.array((types['NO_SFC_ECHO'],types['STRATIFORM'],types['CONVECTIVE'],
                                  types['UNCERTAIN'],types['ISO_CONV_CORE'],types['ISO_CONV_FRINGE'],
                                  types['WEAK_ECHO']))
    rtVar.flag_meanings = np.array(['NO_SFC_ECHO   STRATIFORM   CONVECTIVE   UNCERTAIN   ISO_CONV_CORE   ISO_CONV_FRINGE   WEAK_ECHO'])
    rtVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
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

    rtVar[0,:,:] = rt

    # close file
    ncid.close()

def write_reflCFnetcdf(ncname,refl,
                       deepcoszero,shallowconvmin,minZdiff,truncZconvthres,dBZformaxconvradius,
                       weakechothres,backgrndradius,maxConvRadius,minsize,startslope,maxsize,
                       title,institution,source_uw,references,comment,
                       timeVal,xVal,yVal,latVal,lonVal,gmVal,
                       lat_origin,lon_origin,missing_value):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S")

    # convert timeVal to date and time
    date = dt.datetime.utcfromtimestamp(timeVal[0])
    datetime = date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # open a new netcdf file (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    time = ncid.createDimension('time',None) # None implies UNLIMITED
    y = ncid.createDimension('y',refl.shape[1])
    x = ncid.createDimension('x',refl.shape[0])

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

    reflVar = ncid.createVariable('rain_type',np.int32,('time','y','x'),zlib=True, \
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

    # raintype
    reflVar.units = 'none'
    reflVar.long_name = 'rain_type_classification'
    #reflVar.grid_mapping = 'azimuthal_equidistant'
    reflVar.coordinates = 'lon lat'
    reflVar.grid_mapping = 'grid_mapping'
    reflVar.ancillary_variables = 'rt_Z_th rt_R_bg rt_a rt_b rt_R_conv rt_Z_conv rt_Z_weak rt_Z_shallow rt_A_low rt_A_med rt_A_high'

    # create global attributes
    ncid.Conventions = "CF-1.0"
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    ncid.source_uw = source_uw
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

    reflVar[0,:,:] = refl

    # close file
    ncid.close()

