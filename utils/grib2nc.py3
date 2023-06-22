#!/usr/bin/python3

import os
import sys

if len(sys.argv) != 4:
    print('Usage: {} [inpath] [infile] [outpath]'.format(sys.argv[0]))
    sys.exit()
else:
    inpath = sys.argv[1]
    infile = sys.argv[2]
    outpath = sys.argv[3]

#print('inpath  =',inpath)
#print('infile  =',infile)
#print('outpath =',outpath)

if infile.endswith('grib2'):
    outfile = infile.replace('grib2','nc')
    #print('outfile =',outfile)
    command = 'wgrib2 '+inpath+'/'+infile+' -netcdf '+outpath+'/'+outfile
    #print(command)
    os.system(command)

