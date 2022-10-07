import os
import netCDF4 as nc4
import numpy as np

def ZtoDBZ(z):
  return 10*np.log10(z)

def DBZtoZ(dbz):
  return 10**(0.1*dbz)

indir = '/home/disk/shear2/brodzik/PNNL/po-lun'

for fname in os.listdir(indir):

    #if fname.endswith('nc'):
    if fname.startswith('acme_cmdv_test1.2008-05-2'):

        print fname

        # open fname
        ncid = nc4.Dataset(indir+'/'+fname,'a')
        
        # read var of interest (refl)
        refl = ncid.variables['DBZE_CS'][:]
        (nlev,ncol,nlat,nlon) = refl.shape
        missing = ncid.variables['DBZE_CS'].missing_value
        refl[(refl==missing)] = np.nan

        # linearize refl
        refl_lin = DBZtoZ(refl)

        # create array for averaged data
        refl_lin_avg = np.zeros((nlev,nlat,nlon),dtype=float)

        # average columns for each grid point
        refl_lin_avg = np.nanmean(refl_lin, axis=1)

        # turn linear values back into log values and change nans back to missing
        refl_avg = ZtoDBZ(refl_lin_avg)
        refl_avg[(np.isnan(refl_avg))] = missing

        # create new 3D refl variable
        refl_id = ncid.createVariable('DBZE','f4',('lev','lat','lon'),zlib=True)
        refl_id.units = "dBZe"
        refl_id.long_name = "Column averaged radar dBZe in each gridpoint"
        refl_id.missing_value = missing
        refl_id[:,:,:] = refl_avg

        # close input file
        ncid.close()


