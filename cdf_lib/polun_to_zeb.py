import os
import netCDF4 as nc4
import numpy as np
import math
import datetime as dt

indir = '/home/disk/shear2/brodzik/PNNL/po-lun/20170629'
outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/trmm_sim'

refl_out_missing = -999

MtoKM = 1./1000.

for fname in os.listdir(indir):

    if fname.endswith('nc'):
    #if fname.startswith('acme_cmdv_test1.2008-05'):

        print fname
        
        # get date and time and convert to seconds since Jan 1, 1970
        tmp = fname.split('.')
        date_time = tmp[0]
        [year,month,day,hour] = date_time.replace('-',' ').replace('_',' ').split()
        t = dt.datetime(int(year),int(month),int(day),int(hour))
        basetime = int( (t - dt.datetime(1970,1,1)).total_seconds() )
        
        # open fname
        ncid = nc4.Dataset(indir+'/'+fname,'r')
        
        # read vars of interest
        time_in = ncid.variables['time'][:]
        ntime = time_in.size
        lat_in = ncid.variables['lat'][:]
        nlat = lat_in.size
        lon_in = ncid.variables['lon'][:]
        nlon = lon_in.size
        lev_in = ncid.variables['lev'][:]
        lev_in = lev_in * MtoKM
        nlev = lev_in.size
        
        pres_in = ncid.variables['PS'][:]
        temp_in = ncid.variables['T'][:]

        refl_avg_in = ncid.variables['DBZE_avg'][:]
        missing_avg = ncid.variables['DBZE_avg'].missing_value
        refl_avg_in[refl_avg_in==missing_avg] = refl_out_missing

        refl_max_in = ncid.variables['DBZE_max'][:]
        missing_max = ncid.variables['DBZE_max'].missing_value
        refl_max_in[refl_max_in==missing_max] = refl_out_missing

        counts_in = ncid.variables['valid_counts'][:]
        missing_counts = ncid.variables['valid_counts'].missing_value

        # close input file
        ncid.close()

        # convert lons to values between plus and minus 180
        lon_in = lon_in - 360.

        # create new zeb-formatted netcdf file
        ncid = nc4.Dataset(outdir+'/'+fname,'w',format='NETCDF4')
        
        # create dimensions
        time = ncid.createDimension('time',None) # None implies UNLIMITED
        z = ncid.createDimension('altitude',nlev)
        y = ncid.createDimension('latitude',nlat)
        x = ncid.createDimension('longitude',nlon)

        # create variables
        bt = ncid.createVariable('base_time',np.int32 )
        to = ncid.createVariable('time_offset',np.float32,('time',) )
        alt = ncid.createVariable('altitude',np.float32,('altitude',) )
        lat = ncid.createVariable('latitude',np.float32,('latitude',) )
        lon = ncid.createVariable('longitude',np.float32,('longitude',) )
        pres = ncid.createVariable('sfc_pres',np.float32,('time','latitude','longitude',) )
        temp = ncid.createVariable('temp',np.float32,('time','altitude','latitude','longitude',) )
        refl_avg = ncid.createVariable('refl_avg',np.float32,('time','altitude','latitude','longitude',) )
        refl_max = ncid.createVariable('refl_max',np.float32,('time','altitude','latitude','longitude',) )
        counts = ncid.createVariable('valid_counts',np.float32,('time','altitude','latitude','longitude',) )

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
        refl_avg.units = 'dBZ'
        refl_avg.long_name = 'Simulated average TRMM-PR Reflectivity'
        refl_avg.missing_value = refl_out_missing
        refl_max.units = 'dBZ'
        refl_max.long_name = 'Simulated maximum TRMM-PR Reflectivity'
        refl_max.missing_value = refl_out_missing
        counts.units = 'none'
        counts.long_name = 'Number of valid values (0-50) in each gridpoint'
        counts.missing_value = missing_counts

        # write vars to file
        bt[:] = basetime
        to[:] = 0
        alt[:] = lev_in[:]
        lat[:] = lat_in[:]
        lon[:] = lon_in[:]
        pres[:,:,:] = pres_in
        temp[:,:,:,:] = temp_in
        refl_avg[:,:,:,:] = refl_avg_in
        refl_max[:,:,:,:] = refl_max_in
        counts[:,:,:,:] = counts_in

        # close file
        ncid.close()

        
