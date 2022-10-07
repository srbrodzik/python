import os
import numpy as np
import fnmatch as fnm

region = 'ASIA'
sub_region = 'CHINA'
#limits=[15.0,90.0,45.0,150.0]   # initial run
limits=[15.0,90.0,50.0,150.0]
#indir = '/home/disk/shear2/brodzik/gpm/ERAi_analysis/data'
indir = '/home/disk/shear2/brodzik/IDL/ERAi_analysis/data'
seasons = ['DJF','MAM','JJA','SON']
events = ['DWC','WCC']
thresholds = ['str','mod']
outdir=indir

for season in seasons:
    print season
    if season == 'DJF':
        months = ['12','01','02']
    elif season == 'MAM':
        months = ['03','04','05']
    elif season == 'JJA':
        months = ['06','07','08']
    elif season == 'SON':
        months = ['09','10','11']
    else:
        print 'season='+season+' is invalid . . . exiting'
        exit()
    for event in events:
        print event
        for threshold in thresholds:
            print threshold
            #read input file and save records that are within limits
            for fname in os.listdir(indir):
                if fnm.fnmatch(fname,region+'_'+event+'_'+threshold+'_'+season+'_uw.txt'):
                    print fname
                    f_in = open(indir+'/'+fname,'r')
                    outfile = sub_region+'_'+event+'_'+threshold+'_'+season+'_uw.txt'
                    f_out = open(outdir+'/'+outfile,'w')
                    ncases = f_in.readline()
                    f_out.write(ncases)
                    for line in f_in:
                        #line = line.strip()
                        columns = line.split()
                        lon = float(columns[3])
                        lat = float(columns[4])
                        if lat>=limits[0] and lat<=limits[2] and lon>=limits[1] and lon<=limits[3]:
                            f_out.write(line)
                    f_in.close()
                    f_out.close()
