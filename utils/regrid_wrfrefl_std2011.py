# to execute it interactively
# request nodes: $ salloc -N 1 -p debug  -t 00:30:00
# execute:  $ srun -n ncores python regrid_wrfrefl.py
#    where ncores is the number of cores to request
#
# to execute it in batch
# enter $sh exec_on_cndebug.sh
# enter python code: regrid_wrfrefl.py
# messages  will be written to regrid_wrfrefl.py.out 
from __future__ import print_function
import os.path
import numpy as np
from netCDF4 import Dataset, chartostring
from wrf import getvar, interplevel
from mpi4py import MPI

#outdir = '/scratch3/scratchdirs/feng045/usa/wrfout/wrf.out.run.LGdm.STD2011/'
outdir = '/scratch3/scratchdirs/feng045/usa/wrf/wrf.out.run.LGdm.STD2011/reflectivity/'
outbasename = 'wrfreg_REFL_'

with open('filelist_STD201108.txt') as f:
    lines = f.read().splitlines()

# Define vertical levels
lev = np.arange(500,20001,500)

size = MPI.COMM_WORLD.Get_size()
#print(size)
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
name = MPI.Get_processor_name()
ncores = comm.Get_size()
nlines = len(lines)
print('Number of CPU cores:' + str(ncores))
print('Number of files:' + str(nlines))

# Check modulus (i.e. remainder of nlines / ncores)
nfremain = float(nlines) % float(ncores)
# If remainder > 0 (i.e., non-divisible by ncores), increase number of lines by remainder
# such that nlines are always divisible by ncores
if (nfremain > 0): 
	nlines = nlines + ncores - nfremain

nfpercore = long(nlines / ncores)
print('Number of files per core:' + str(nfpercore))

# Send a batch of files to each rank (CPU #)
for ifile in range(rank*nfpercore,(rank+1)*nfpercore):

	# filename
	filein = lines[ifile]
	# Get date/time string from filename
	# assuming wrfout file has this standard format: wrfout_d0x_yyyy-mm-dd_hh:mm:ss
	fname = os.path.basename(filein)
	ftime = fname[11:]
	fileout = outdir + outbasename + ftime + '.nc'
	print(fileout)
	print(lines[rank])

	# Read WRF file
	ncfile = Dataset(lines[rank])
	#times = getvar(ncfile, "Times")
	times = ncfile.variables['Times'][:]
	Timesstr = np.array2string(np.squeeze(chartostring(times)))
	print(Timesstr)
	strlen = len(Timesstr[1:-1])
	DX = getattr(ncfile, 'DX')
	DY = getattr(ncfile, 'DY')

	z = getvar(ncfile, "z")
	XLAT = getvar(ncfile, "XLAT")
	XLONG = getvar(ncfile, "XLONG")
	lat1d = XLAT[:,0]
	lon1d = XLONG[0,:]
	dims = np.shape(XLAT)
	ny = dims[0]
	nx = dims[1]
	#p = getvar(ncfile, "pressure")
	refl = getvar(ncfile, "REFL_10CM")

	# Convert reflectivity to linear unit
	refllin = 10.0**(refl/10.0)
        
	# Interpolate to fixed vertical level
	refllinreg = np.zeros((len(lev),ny,nx), dtype=np.float)
	for i in range(0,len(lev)):
		refllinreg[i,:,:] = interplevel(refllin, z, lev[i])

	# Convert reflectivity to log (dBZ) unit
	reflreg = 10.0 * np.log10(refllinreg)

	# Replace bad values with NAN
	fillval = np.nan
	reflreg[np.where(reflreg > 100)] = fillval
        

	# Write output to netCDF file
	regdata = Dataset(fileout, 'w', format='NETCDF4_CLASSIC')
	regdata.createDimension('time', None)
	regdata.createDimension('height', len(lev))
	regdata.createDimension('lat', ny)
	regdata.createDimension('lon', nx)
	regdata.createDimension('DateStrLen',strlen)

	# Set global attributes
	regdata.setncattr('DX', DX)
	regdata.setncattr('DY', DY)
	
	# Note: variable object (e.g. lat) must not have a variable named the same
	#       otherwise the output would be empty.
	Times = regdata.createVariable('Times', 'S1', ('time','DateStrLen'))

#	lon = regdata.createVariable('lon', 'f4', 'lon', zlib=True)
#	lon.units = 'degrees_east'
#	lon.long_name = 'longitude'
#	
#	lat = regdata.createVariable('lat', 'f4', 'lat', zlib=True)
#	lat.units = 'degrees_north'
#	lat.long_name = 'latitude'

        lon2d = regdata.createVariable('lon2d', 'f4', ('lat','lon'), zlib=True, complevel=5)
        lon2d.units = 'degrees_east'
        lon2d.long_name = 'longitude'

        lat2d = regdata.createVariable('lat2d', 'f4', ('lat','lon'), zlib=True, complevel=5)
        lat2d.units = 'degrees_north'
        lat2d.long_name = 'latitude'

	height = regdata.createVariable('height', 'f4', 'height', zlib=True)
	height.units = 'm'
	height.long_name = 'height level'
	
	REFL = regdata.createVariable('REFL', 'f4', ('time','height','lat','lon'), zlib=True, complevel=5, fill_value=fillval)
	REFL.units = 'DBZ'
	REFL.description = 'REFL_10CM'
	
	Times[0,:] = times
#	lon[:] = lon1d[:]
#	lat[:] = lat1d[:]
	lon2d[:,:] = XLONG[:,:]
	lat2d[:,:] = XLAT[:,:]
	height[:] = lev[:]
	REFL[0,:,:,:] = reflreg[:,:,:]
	
	regdata.close()

