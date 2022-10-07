import os
import netCDF4 as nc4
import numpy as np

def ZtoDBZ(z):
  return 10*np.log10(z)

def DBZtoZ(dbz):
  return 10**(0.1*dbz)

indir = '/home/disk/shear2/brodzik/PNNL/po-lun/20170629'
missing_int = -99

for fname in os.listdir(indir):

    if fname.endswith('nc'):
    #if fname.startswith('acme_cmdv_test1.2008-05-2'):

        print fname

        # open fname
        ncid = nc4.Dataset(indir+'/'+fname,'a')
        
        # read var of interest (refl)
        #refl = ncid.variables['DBZE_CS'][:]
        refl = np.array(ncid.variables['DBZE_CS'])
        (ntime,nlev,ncol,nlat,nlon) = refl.shape
        missing = ncid.variables['DBZE_CS'].missing_value
        refl[(refl==missing)] = np.nan

        # linearize refl
        refl_lin = DBZtoZ(refl)

        # create arrays for outputs
        refl_lin_avg = np.zeros((ntime,nlev,nlat,nlon),dtype=float)
        refl_lin_max = np.zeros((ntime,nlev,nlat,nlon),dtype=float)
        valid_counts = np.zeros((ntime,nlev,nlat,nlon),dtype=int)

        # average columns for each grid point
        refl_lin_avg = np.nanmean(refl_lin, axis=2)

        # find max of columns for each grid point
        refl_lin_max = np.nanmax(refl_lin, axis=2)

        # count number of non-nan values at each gridpoint
        #valid_counts = np.count_nonzero(~np.isnan(refl),axis=2)
        for itime in range(0,ntime):
          for ilev in range(0,nlev):
            for ilat in range(0,nlat):
              for ilon in range(0,nlon):
                valid_counts[itime,ilev,ilat,ilon] = np.count_nonzero(~np.isnan(refl[itime,ilev,:,ilat,ilon]))

        # turn linear values back into log values and change nans back to missing
        refl_avg = ZtoDBZ(refl_lin_avg)
        refl_avg[(np.isnan(refl_avg))] = missing
        refl_max = ZtoDBZ(refl_lin_max)
        refl_max[(np.isnan(refl_max))] = missing

        # create new 3D variables
        refl_avg_id = ncid.createVariable('DBZE_avg','f4',('time','lev','lat','lon'),zlib=True)
        refl_avg_id.units = "dBZe"
        refl_avg_id.long_name = "Column averaged radar dBZe in each gridpoint"
        refl_avg_id.missing_value = missing
        refl_avg_id[:,:,:,:] = refl_avg
        
        refl_max_id = ncid.createVariable('DBZE_max','f4',('time','lev','lat','lon'),zlib=True)
        refl_max_id.units = "dBZe"
        refl_max_id.long_name = "Column max of radar dBZe in each gridpoint"
        refl_max_id.missing_value = missing
        refl_max_id[:,:,:,:] = refl_max

        valid_counts_id = ncid.createVariable('valid_counts',np.int32,('time','lev','lat','lon'),zlib=True)
        valid_counts_id.units = "none"
        valid_counts_id.long_name = "Number of valid values (0-50) in each gridpoint"
        valid_counts_id.missing_value = missing_int
        valid_counts_id[:,:,:,:] = valid_counts

        # close input file
        ncid.close()


        
