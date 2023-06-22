# Brody Fuchs, CSU, Sept 27 2013
# brfuchs@atmos.colostate.edu

# This code is designed to take in a netcdf file and "delete" variables
# this done by reading in the netcdf file then writing back everything
# to the new file except for the unwanted variables
# If someone else is running this script then will need to change 
# filenames and things


import numpy as np
from netCDF4 import Dataset
import datetime as dt
import glob
import os

del_var = ['lma_'] # string with variable(s) to be deleted
del_dim = ['s', 'f'] # associated dimensions to be deleted

dir = '/data2/brfuchs/research'
region = 'DC'
input = open('%s/tracking/processdate.txt'%dir).readlines()

zipold = 1

year,month,date,n_days = np.int(input[0]),np.int(input[1]),np.int(input[2]),np.int(input[3])

in_dir = '/data2/brfuchs/research/radartest/aero_nc' # files that are getting taken in
out_dir = '/data2/brfuchs/research/radartest/nldn_nc' # dir that the files go out to
dlist = sorted(glob.glob('%s/%d/%d????'%(in_dir,year,year)))
#dlist.sort()
sdlist = [os.path.basename(i) for i in dlist] # list of directories with date names

april = dt.datetime(year,4,1)
a_dir = dt.datetime.timetuple(april)[7]

starttime = dt.datetime(year, month, date)
s_dir = dt.datetime.timetuple(starttime)[7]
endtime = starttime+dt.timedelta(days=n_days-1)
e_dir = dt.datetime.timetuple(endtime)[7]
s = s_dir-a_dir
e = e_dir-a_dir+1


for iday in range(s,e):
    print 'day is %s'%sdlist[iday]
    files = sorted(glob.glob('%s/mosaic3d_%s_*.netcdf'%(dlist[iday],region)))
    for file in files:
        try:
            grp = Dataset(file)

            var_keys = grp.variables.keys()
            dim_keys = grp.dimensions.keys()
            att_keys = grp.ncattrs()
            # check for variables that will survive

	# GOOD VARS NEEDS TO GET FIXED!!!!!!!!!

            good_vars = [i for i in var_keys if i.find(del_var[0]) == -1] # if string del_var isnt found then keep it
            good_dims = [i for i in dim_keys if i not in del_dim] # getting rid of del_dim from dims

	    if len(del_var) > 1:
	        for k in del_var:
                    good_vars = [i for i in good_vars if i.find(k) == -1] # if string del_var isnt found then keep it


            outvar = {} # actual variable
            outatt = {} # attributes
            outdim = {} # dimensions
            outvardim = {} # variable dimensions
            outvartype = {} # variable type
            for var in good_vars:
                outvar[var] = grp.variables[var][:]
                outvardim[var] = grp.variables[var].dimensions
                outvartype[var] = grp.variables[var].dtype
            for att in att_keys:
                outatt[att] = grp.getncattr(att)
            for dim in good_dims:
                outdim[dim] = len(grp.dimensions[dim])

            grp.close()
        except RuntimeError:
            pass

        # GOT ALL THE INFO WE NEED, NOW MAKE NEW FILE
	f_base = os.path.basename(file)
        outfile = '%s/%d/%s/%s'%(out_dir, year, sdlist[iday],f_base)
        w_grp = Dataset(outfile, 'w')
        # FILE OPENED, NOW WRITE THE DIMENSIONS
        for i in outdim.keys():
            w_grp.createDimension(i, outdim[i])
        for i in range(len(outvar.keys())):
            dummy = w_grp.createVariable(outvar.keys()[i], outvartype[outvartype.keys()[i]], outvardim[outvardim.keys()[i]])  
            dummy[:] = outvar[outvar.keys()[i]]
        for i in outatt.keys():
            w_grp.setncattr(i, outatt[i])


        w_grp.close()

    if zipold == 1:
	print 'zipping old %s files'%sdlist[iday]
	os.system('gzip %s/%d/%s/*%s*.netcdf'%(in_dir,year,sdlist[iday],region))


