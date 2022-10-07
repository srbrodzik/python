from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import netCDF4 as nc4
import os
import logging as log

# input and output directories
fileDirIn = '/home/disk/mjo/dynamo/data.server/zebra/QCed/spolka/sur_1km_legacy/20111016/'

# modify standard_name for these vars
modVar = ['REFL','VEL']

# new standard_names for modVar
modVal = ['equivalent_reflectivity_factor',  \
          'radial_velocity_of_scatterers_away_from_instrument']

# remove standard_name attribute for these var
delVar = ['WIDTH','ZDR','LDR','RHOHV','KDP','PID']

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for fname in os.listdir(fileDirIn):
    if fname.endswith('nc'):
        log.info( "file = {}".format(fname) )

        # Open input file
        ncid = nc4.Dataset(str(fileDirIn + fname), 'a')

        # Modify standard_name values
        for var in modVar:
            index = modVar.index(var)
            log.info( "var = {} standard_name = {}".format(var,modVal[index]) )
            ncid.variables[var].standard_name = modVal[index]

        # Remove standard_name attributes
        for var in delVar:
            del ncid.variables[var].standard_name
            log.info( "var = {} standard_name removed".format(var) )

        # Close netcdf file
        ncid.close()
