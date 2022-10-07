import netCDF4 as nc4
#import os
import numpy as np
import time as tm
import datetime as dt
import calendar

def readTRMMascii(dataDir,prefix,levels,delta_z,missingVal):

    # get orbit and time from prefix
    orbit = prefix[5:10]
    date = prefix[11:19]
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    time = prefix[20:26]
    hours = time[0:2]
    mins = time[2:4]
    secs = time[4:6]

    # convert  date and time into seconds since 1/1/1970 at 00:00
    timeStr = year + ' ' + month + ' ' + day + ' ' + hours + ' ' + mins + ' ' + secs
    baseTime = calendar.timegm(tm.strptime(timeStr,'%Y %m %d %H %M %S') )
    
    # change working directory
    #os.chdir(dataDir)

    # read number of columns (lons)
    ncol = 0
    f = open(prefix+'gridColumns.txt','r')
    ncol = int(f.read())
    f.close()
    
    # read number of rows (lats)
    nrow = 0
    f = open(prefix+'gridRows.txt','r')
    nrow = int(f.read())
    f.close()

    # read lats and lons
    lats1D = np.loadtxt(prefix+'gridLatitude.txt',dtype=float)
    lats3D = np.reshape( lats1D, (int(ncol),int(nrow)) )
    lons1D = np.loadtxt(prefix+'gridLongitude.txt',dtype=float)
    lons3D = np.reshape( lons1D, (int(ncol),int(nrow)) )

    # read raintype
    rt1D = np.loadtxt(prefix+'gridRainType.txt',dtype=float)
    rt1D[rt1D == -88] = missingVal
    rt1D[np.isnan(rt1D)] = missingVal
    raintype = np.reshape( rt1D, (int(ncol),int(nrow)) )
    
    # read swath
    swath1D = np.loadtxt(prefix+'noSwath.txt',dtype=float)
    swath1D[swath1D <= missingVal] = missingVal
    swath = np.reshape( swath1D, (int(ncol),int(nrow)) )

    # read reflectivity
    refl1D = np.loadtxt(prefix+'gridReflectivity3D3L.txt',dtype=float)
    refl1D[refl1D <= missingVal] = missingVal
    refl = np.reshape( refl1D, (levels,int(ncol),int(nrow)) )

    # transpose arrays as needed
    n_lats = nrow
    lats = lats3D[0,:]
    n_lons = ncol
    lons = lons3D[:,0]
    raintype = raintype.transpose((1,0))
    swath = swath.transpose((1,0))
    refl = refl.transpose((0,2,1))
    alts = np.linspace(0,levels-1,num=levels) * delta_z

    return (baseTime,lats,lons,alts,raintype,swath,refl)

def readGPMascii(dataDir,prefix,levels,delta_z,missingVal):

    # get orbit and time from prefix
    orbit = prefix[8:14]
    date = prefix[15:23]
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    time = prefix[24:30]
    hours = time[0:2]
    mins = time[2:4]
    secs = time[4:6]
    # missingVal=-99  #from idl code
    missingValInt=-9999
    missingValByte=255
    noRainValInt=-1111

    # convert  date and time into seconds since 1/1/1970 at 00:00
    timeStr = year + ' ' + month + ' ' + day + ' ' + hours + ' ' + mins + ' ' + secs
    baseTime = calendar.timegm(tm.strptime(timeStr,'%Y %m %d %H %M %S') )
    
    # change working directory
    #os.chdir(dataDir)

    # read number of columns (lons)
    ncol = 0
    f = open(prefix+'gridColumns.txt','r')
    ncol = int(f.read())
    f.close()
    
    # read number of rows (lats)
    nrow = 0
    f = open(prefix+'gridRows.txt','r')
    nrow = int(f.read())
    f.close()

    # USING np.loadtxt INSTEAD OF np.genfromtxt BECAUSE IT'S FASTER
    # read lats and lons
    lats1D = np.loadtxt(prefix+'gridLatitude.txt',dtype=float)
    lats3D = np.reshape( lats1D, (int(ncol),int(nrow)) )
    lons1D = np.loadtxt(prefix+'gridLongitude.txt',dtype=float)
    lons3D = np.reshape( lons1D, (int(ncol),int(nrow)) )

    # read raintype
    rt1D = np.loadtxt(prefix+'gridRainType.txt',dtype=float)
    rt1D[rt1D == -1111] = missingVal
    rt1D[np.isnan(rt1D)] = missingVal
    rainType = np.reshape( rt1D, (int(ncol),int(nrow)) )
    
    # read raw raintype
    rtr1D = np.loadtxt(prefix+'gridRainTypeORG.txt',dtype=float)
    rtr1D[np.isnan(rtr1D)] = missingVal
    rainTypeRaw = np.reshape( rtr1D, (int(ncol),int(nrow)) )
    
    # read phase type
    pt1D = np.loadtxt(prefix+'gridNearSurfPhase.txt',dtype=float)
    pt1D[pt1D < 0] = missingVal    # missing vals set to -9999
    pt1D[pt1D == missingValByte] = missingVal
    pt1D[np.isnan(pt1D)] = missingVal
    phaseType = np.reshape( pt1D, (int(ncol),int(nrow)) )
    
    # read raw phase type
    ptr1D = np.loadtxt(prefix+'gridNearSurfPhaseORG.txt',dtype=float)
    ptr1D[ptr1D < 0] = missingVal
    ptr1D[ptr1D == missingValByte] = missingVal
    ptr1D[np.isnan(ptr1D)] = missingVal
    phaseTypeRaw = np.reshape( ptr1D, (int(ncol),int(nrow)) )
    
    # read shallow rain type
    shrn1D = np.loadtxt(prefix+'gridShallowRainType.txt',dtype=float)
    shrn1D[shrn1D < 0] = missingVal    # missing vals = -9999 & no rain val = -1111
    shallowRainType = np.reshape( shrn1D, (int(ncol),int(nrow)) )
    
    # read near surface rain
    nsr1D = np.loadtxt(prefix+'gridNearSurfRain.txt',dtype=float)
    nsr1D[nsr1D < 0] = missingVal
    nearSurfRain = np.reshape( nsr1D, (int(ncol),int(nrow)) )
    
    # read bright band width
    bbw1D = np.loadtxt(prefix+'gridWidthBB.txt',dtype=float)
    bbw1D[bbw1D < 0] = missingVal    # missing vals = -9999. & no rain val = -1111.1
    bbWidth = np.reshape( bbw1D, (int(ncol),int(nrow)) )
    
    # read bright band height
    bbh1D = np.loadtxt(prefix+'gridHeightBB.txt',dtype=float)
    bbh1D[bbh1D < 0] = missingVal    # missing vals = -9999 & no rain val = -1111.1
    bbHeight = np.reshape( bbh1D, (int(ncol),int(nrow)) )
    
    # read swath
    swath1D = np.loadtxt(prefix+'noSwath.txt',dtype=float)
    swath1D[swath1D <= missingVal] = missingVal
    swath = np.reshape( swath1D, (int(ncol),int(nrow)) )

    # read reflectivity
    refl1D = np.loadtxt(prefix+'gridReflectivity3D3L.txt',dtype=float)
    refl1D[refl1D <= missingVal] = missingVal
    refl = np.reshape( refl1D, (levels,int(ncol),int(nrow)) )

    # read hdf file name
    f = open(prefix+'hdfFileName.txt','r')
    hdfFileName = f.read()
    f.close()
    hdfFileName = hdfFileName+'.HDF5'
             
    # transpose arrays as needed
    n_lats = nrow
    lats = lats3D[0,:]
    n_lons = ncol
    lons = lons3D[:,0]
    rainType = rainType.transpose((1,0))
    rainTypeRaw = rainTypeRaw.transpose((1,0))
    phaseType = phaseType.transpose((1,0))
    phaseTypeRaw = phaseTypeRaw.transpose((1,0))
    shallowRainType = shallowRainType.transpose((1,0))
    nearSurfRain = nearSurfRain.transpose((1,0))
    bbWidth = bbWidth.transpose((1,0))
    bbHeight = bbHeight.transpose((1,0))
    swath = swath.transpose((1,0))
    refl = refl.transpose((0,2,1))
    alts = np.linspace(0,levels-1,num=levels) * delta_z

    return (baseTime,lats,lons,alts,rainType,rainTypeRaw,phaseType,phaseTypeRaw,
            shallowRainType,nearSurfRain,bbWidth,bbHeight,swath,refl,hdfFileName)

def writeZeb_TRMMnetcdf(dataDir,prefix,baseTime,lats,lons,alts,raintype,swath,refl,
                    title,institution,missingVal):

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncname = prefix + 'nc'
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    timeDim = ncid.createDimension('time',None) # None implies UNLIMITED
    altDim = ncid.createDimension('altitude',refl.shape[0])
    latDim = ncid.createDimension('latitude',refl.shape[1])
    lonDim = ncid.createDimension('longitude',refl.shape[2])

    # create variables
    bt = ncid.createVariable('base_time',np.float64 )
    to = ncid.createVariable('time_offset',np.float32,('time',) )
    la = ncid.createVariable('latitude',np.float32,('latitude') )
    lo = ncid.createVariable('longitude',np.float32,('longitude',) )
    al = ncid.createVariable('altitude',np.float32,('altitude',) )
    rt = ncid.createVariable('rain_type',np.float32,('time','latitude','longitude',), \
                                   fill_value=missingVal )
    sw = ncid.createVariable('swath',np.float32,('time','latitude','longitude',), \
                                   fill_value=missingVal )
    re = ncid.createVariable('refl',np.float32,('time','altitude','latitude','longitude',), \
                                   fill_value=missingVal )

    # create variable attributes
    bt.units = 'seconds since 1970-01-01 00:00:00 +0000'
    to.units = 'seconds since base_time'
    la.units = 'degrees_north'
    lo.units = 'degrees_east'
    al.units = 'km'
    rt.units = 'none'
    rt.long_name = 'Rain Type'
    rt.stratiform = '10'
    rt.convective = '20'
    rt.other = '30'
    sw.units = 'none'
    sw.long_name = 'GPM-Ku coverage area'
    re.units = 'dBZ'
    re.long_name = 'GPM_Ku Reflectivity'

    # create global attributes
    ncid.title = title
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    #ncid.source = source
    #ncid.references = references
    #ncid.comment = comment

    # write vars to file
    bt[:] = baseTime
    to[:] = 0
    la[:] = lats
    lo[:] = lons
    al[:] = alts
    rt[0,:,:] = raintype
    sw[0,:,:] = swath
    re[0,:,:,:] = refl

    #close file
    ncid.close()

def writeZeb_GPMnetcdf(dataDir,prefix,baseTime,lats,lons,alts,rainType,rainTypeRaw,
                       phaseType,phaseTypeRaw,shallowRainType,nearSurfRain,bbWidth,
                       bbHeight,swath,refl,hdfFileName,title,institution,missingVal):

    missingValInt=-9999
    missingValByte=255
    noRainValInt=-1111

    # get current time
    currentTime = tm.strftime("%m/%d/%Y %H:%M:%S");

    # open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
    ncname = prefix + 'nc'
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    # create dimensions
    timeDim = ncid.createDimension('time',None) # None implies UNLIMITED
    altDim = ncid.createDimension('altitude',refl.shape[0])
    latDim = ncid.createDimension('latitude',refl.shape[1])
    lonDim = ncid.createDimension('longitude',refl.shape[2])

    # create variables
    bt = ncid.createVariable('base_time',np.float64 )
    to = ncid.createVariable('time_offset',np.float32,('time',) )
    la = ncid.createVariable('latitude',np.float32,('latitude') )
    lo = ncid.createVariable('longitude',np.float32,('longitude',) )
    al = ncid.createVariable('altitude',np.float32,('altitude',) )
    rt = ncid.createVariable('rain_type',np.float32,('time','latitude','longitude',), \
                             fill_value=missingVal )
    rtr = ncid.createVariable('rain_type_raw',np.float32,('time','latitude','longitude',), \
                              fill_value=missingVal )
    pt = ncid.createVariable('phase_type',np.float32,('time','latitude','longitude',), \
                             fill_value=missingVal )
    ptr = ncid.createVariable('phase_type_raw',np.float32,('time','latitude','longitude',), \
                              fill_value=missingVal )
    shrn = ncid.createVariable('shallow_rain_type',np.float32,('time','latitude','longitude',), \
                               fill_value=missingVal )
    nsr = ncid.createVariable('near_surf_rain',np.float32,('time','latitude','longitude',), \
                              fill_value=missingVal )
    bbw = ncid.createVariable('width_bb',np.float32,('time','latitude','longitude',), \
                              fill_value=missingVal )
    bbh = ncid.createVariable('height_bb',np.float32,('time','latitude','longitude',), \
                              fill_value=missingVal )
    sw = ncid.createVariable('swath',np.float32,('time','latitude','longitude',), \
                             fill_value=missingVal )
    re = ncid.createVariable('refl',np.float32,('time','altitude','latitude','longitude',), \
                             fill_value=missingVal )

    # create variable attributes
    bt.units = 'seconds since 1970-01-01 00:00:00 +0000'
    to.units = 'seconds since base_time'
    la.units = 'degrees_north'
    lo.units = 'degrees_east'
    al.units = 'km'

    rt.units = 'none'
    rt.long_name = 'Rain Type'
    rt.stratiform = 1
    rt.convective = 2
    rt.other = 3
    rt.missing_value = missingVal
             
    rtr.units = 'none'
    rtr.long_name = 'Rain Type Raw'
    rtr.stratiform_range = '[10000000,19999999]'
    rtr.convective_range = '[20000000,29999999]'
    rtr.other_range = '[30000000,39999999]'
    rtr.no_rain_value = noRainValInt
    rtr.missing_value = missingValInt
             
    pt.units = 'none'
    pt.long_name = 'Near Surface Phase'
    pt.solid = 0
    pt.mixed = 1
    pt.liquid = 2
    pt.missing_value = missingVal

    ptr.units = 'none'
    ptr.long_name = 'Near Surface Phase Raw'
    ptr.phaseType_lt_100 = 'tempC = phaseType - 100'
    ptr.phaseType_gt_200 = 'tempC = phaseType - 200'
    ptr.phaseType_eq_100 = 'top of bright band'
    ptr.phaseType_eq_125 = 'between top and peak of BB'
    ptr.phaseType_eq_175 = 'between peak and bottom of BB'
    ptr.phaseType_eq_200 = 'bottom of bright band'
    ptr.missing_value = missingValByte

    shrn.units = 'none'
    shrn.long_name = 'Shallow Rain Flag'
    shrn.no_shallow_rain = 0
    shrn.shallow_isolated_maybe = 10 
    shrn.shallow_isolated_certain = 11
    shrn.shallow_nonisolated_maybe = 20
    shrn.shallow_nonisolated_certain = 21
    shrn.no_rain_value = noRainValInt
    shrn.missing_value = missingVal

    nsr.units = 'mm/hr'
    nsr.long_name = 'Near Surface Rain'
    nsr.missing_value = missingVal

    bbw.units = 'meters'
    bbw.long_name = 'Width of Bright Band'
    bbw.missing_value = missingVal

    bbh.units = 'meters'
    bbh.long_name = 'Height of Bright Band'
    bbh.missing_value = missingVal

    sw.units = 'none'
    sw.long_name = 'GPM-Ku coverage area'
    sw.swath = 0
    sw.no_swath = 1 
    sw.missing_value = missingVal
             
    re.units = 'dBZ'
    re.long_name = 'GPM_Ku Reflectivity'
    re.missing_value = missingVal

    # create global attributes
    ncid.title = title
    ncid.source = hdfFileName
    ncid.institution = institution
    ncid.history = 'File created ' + currentTime
    #ncid.references = references
    #ncid.comment = comment

    # write vars to file
    bt[:] = baseTime
    to[:] = 0
    la[:] = lats
    lo[:] = lons
    al[:] = alts
    rt[0,:,:] = rainType
    rtr[0,:,:] = rainTypeRaw
    pt[0,:,:] = phaseType
    ptr[0,:,:] = phaseTypeRaw
    shrn[0,:,:] = shallowRainType
    nsr[0,:,:] = nearSurfRain
    bbw[0,:,:] = bbWidth
    bbh[0,:,:] = bbHeight
    sw[0,:,:] = swath
    re[0,:,:,:] = refl

    #close file
    ncid.close()
