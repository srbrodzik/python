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
    lats1D = np.genfromtxt(prefix+'gridLatitude.txt',dtype=float)
    lats3D = np.reshape( lats1D, (int(ncol),int(nrow)) )
    lons1D = np.genfromtxt(prefix+'gridLongitude.txt',dtype=float)
    lons3D = np.reshape( lons1D, (int(ncol),int(nrow)) )

    # read raintype
    rt1D = np.genfromtxt(prefix+'gridRainType.txt',dtype=float)
    rt1D[rt1D == -88] = missingVal
    rt1D[np.isnan(rt1D)] = missingVal
    raintype = np.reshape( rt1D, (int(ncol),int(nrow)) )
    
    # read swath
    swath1D = np.genfromtxt(prefix+'noSwath.txt',dtype=float)
    swath1D[swath1D <= missingVal] = missingVal
    swath = np.reshape( swath1D, (int(ncol),int(nrow)) )

    # read reflectivity
    refl1D = np.genfromtxt(prefix+'gridReflectivity3D3L.txt',dtype=float)
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
    alts = np.linspace(0,levels-1,num=80) * delta_z

    return (baseTime,lats,lons,alts,raintype,swath,refl)

def writeZeb_netcdf(dataDir,prefix,baseTime,lats,lons,alts,raintype,swath,refl,
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
