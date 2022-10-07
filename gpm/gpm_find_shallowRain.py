import os
import netCDF4 as nc4
import numpy as np

# User inputs
indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2015/06'
outfile = indir+'/filesWithShallowRain.txt'

# Open output file
fout = open(outfile, 'w')

# Go through all netcdf files in indir
for fname in os.listdir(indir):
    
    if fname.endswith('nc'):

        print fname

        # Read shallow_rain_type field
        ncid = nc4.Dataset(indir+'/'+fname,'r')
        srt = ncid.variables['shallow_rain_type'][:]
        ncid.close()

        # Print out files with shallow rain
        vals = np.unique(srt)
        if np.amax(srt) > 0:
            fout.write("fname = {} and vals = {}\n".format(fname,vals))

# Close output file
fout.close()
