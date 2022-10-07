import os
import netCDF4 as nc4
import numpy as np
import math
import datetime as dt

# function to convert pressure (Pa) to height (meters)
def PtoH(P):
  a = 101325
  b = 2.25577e-05
  c = 5.25588
  return (1.0 - math.pow(P/a,1.0/c))/b

# flip an array in a given dimension
def reverse(a, axis=0): 
  idx = [slice(None)]*len(a.shape)
  idx[axis] = slice(None, None, -1)
  return a[idx]

indir = '/home/disk/shear2/brodzik/PNNL/po-lun'
#outdir = '/home/disk/shear2/brodzik/PNNL/po-lun/test_out'
outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/trmm_sim_z'

refl_out_missing = -999

for fname in os.listdir(indir):

    #if fname.endswith('nc'):
    if fname.startswith('acme_cmdv_test1.2008-05'):

        print fname
        
        # get date and time and convert to seconds since Jan 1, 1970
        tmp = fname.split('.')
        date_time = tmp[1]
        [year,month,day,hour] = date_time.split('-')
        t = dt.datetime(int(year),int(month),int(day),int(hour))
        basetime = int( (t - dt.datetime(1970,1,1)).total_seconds() )
        
        # open fname
        ncid = nc4.Dataset(indir+'/'+fname,'r')
        
        # read vars of interest (refl)
        lat_in = ncid.variables['lat'][:]
        nlat = lat_in.size
        lon_in = ncid.variables['lon'][:]
        nlon = lon_in.size
        lev = ncid.variables['lev'][:]
        nlev = lev.size
        
        pres_in = ncid.variables['PS'][:]
        temp_in = ncid.variables['T'][:]

        refl_in = ncid.variables['DBZE'][:]
        missing = ncid.variables['DBZE'].missing_value
        refl_in[refl_in==missing] = refl_out_missing

        # close input file
        ncid.close()

        # convert lons to values between plus and minus 180
        lon_in = lon_in - 360.
        
        # convert pressure levels to standard heights (in km)
        ht = np.zeros(lev.shape)
        for i in range(0,nlev):
          ht[nlev-i-1] = PtoH(lev[i]*100)/1000.

        # create new zeb-formatted netcdf file
        ncid = nc4.Dataset(outdir+'/'+fname,'w',format='NETCDF4')
        
        # create dimensions
        time = ncid.createDimension('time',None) # None implies UNLIMITED
        z = ncid.createDimension('altitude',lev.size)
        y = ncid.createDimension('latitude',lat_in.size)
        x = ncid.createDimension('longitude',lon_in.size)

        # create variables
        bt = ncid.createVariable('base_time',np.int32 )
        to = ncid.createVariable('time_offset',np.float32,('time',) )
        alt = ncid.createVariable('altitude',np.float32,('altitude',) )
        lat = ncid.createVariable('latitude',np.float32,('latitude',) )
        lon = ncid.createVariable('longitude',np.float32,('longitude',) )
        pres = ncid.createVariable('sfc_pres',np.float32,('time','latitude','longitude',) )
        temp = ncid.createVariable('temp',np.float32,('time','altitude','latitude','longitude',) )
        #refl = ncid.createVariable('refl',np.float32,('time','altitude','latitude','longitude',),fill_value=missing )
        refl = ncid.createVariable('refl',np.float32,('time','altitude','latitude','longitude',) )

        # create variable attributes
        bt.units = 'seconds since 1970-01-01 00:00:00 +0000'
        to.units = 'seconds since base_time'
        alt.units = 'km'
        alt.long_name = 'altitude MSL'
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        pres.units = 'Pa'
        pres.long_name = 'Surface pressure'
        temp.units = 'degK'
        temp.long_name = 'Temperature'
        refl.units = 'dBZ'
        refl.long_name = 'Simulated TRMM-PR Reflectivity'
        refl.missing_value = refl_out_missing

        # write vars to file
        bt[:] = basetime
        to[:] = 0
        alt[:] = ht[:]
        lat[:] = lat_in[:]
        lon[:] = lon_in[:]
        pres[0,:,:] = pres_in[:,:]
        temp[0,:,:,:] = reverse(temp_in,axis=0)
        #refl[0,:,:,:] = refl_in[:,:,:]
        refl[0,:,:,:] = reverse(refl_in,axis=0)

        # close file
        ncid.close()

        
