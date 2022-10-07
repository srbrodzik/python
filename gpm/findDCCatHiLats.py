# Use this code to find extreme events at high latitudes whose centroids are not in the 60-65N
# range but which extend into those latitudes

import os
import glob
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma

region = 'asia'
region_in_caps = region.upper()
event = 'dcc'
event_in_caps = event.upper()
flistDir = '/home/disk/bob/gpm/'+region+'_ku/classify/class_data_v05_uw_hilat/monthly_class_v10s'
flist  = 'ExtremeEvents_'+region_in_caps+'_HiLat_check_'+event_in_caps+'.csv'
cdfdir = '/home/disk/bob/gpm/'+region+'_ku/classify/ex_data_v05_hilat'
maskName = event+'_mask_str'
LAT_THRESH = 60.0
missing = -99
fout = open(flistDir+'/ExtremeEvents_'+region_in_caps+'_HiLat_addto_'+event_in_caps+'.txt','w')

# FOR TESTING
# fname = ''/home/disk/bob/gpm/aka_ku/classify/ex_data_v05_hilat/2014/07/GPM2Ku5_uw3_20140703.203924_to_20140703.204254_001963_AKA.nc'

# Check file associated with each line in flist
with open(flistDir+'/'+flist) as f:
    for line in f:
        line = line.strip('\n')

        # FOR TESTING AKA
        # line = '001963 20140703 203924 1'    Not in hiLat region
        # line = '008157 20150805 225222 1'    In hiLat region
        
        #fout.write('%s\n' % (line))
        print('%s' % (line))

        parts = line.split(' ')
        orbit = parts[0]
        date  = parts[1]
        time  = parts[2]
        maskNumber = int(parts[3])

        year = date[0:4]
        month = date[4:6]

        # find fname
        fname1 = 'GPM2Ku5_uw3_'+date+'.'+time+'_to_*_'+orbit+'_'+region_in_caps+'.nc'
        fname2 = glob.glob(cdfdir+'/'+year+'/'+month+'/'+fname1)
        fname = fname2[0]
        #fout.write('%s\n' % (fname))

        # open fname and get lat, lon and mask
        ncid = nc4.Dataset(fname,'r')
        lat = np.array(ncid.variables['lat'])
        lon = np.array(ncid.variables['lon'])
        mask = np.array(ncid.variables[maskName])
        mask = np.squeeze(mask)
        ncid.close()

        # FOR TESTING
        # unique,counts = np.unique(mask, return_counts=True)
        # dict(zip(unique,counts))
        
        # first check to see if maxLat > LAT_THRESH
        if np.amax(lat) > LAT_THRESH:
            #fout.write('Passes first test\n')
            print('Passes first test\n')

            # then check to see if any part of mask lies >= minLat
            lon2d, lat2d = np.meshgrid(lon,lat)

            lat2d[mask != maskNumber] = missing
            
            if np.amax(lat2d) > LAT_THRESH:
                # event goes into hiLat region -- save filename
                #fout.write('Passes second test\n')
                print('Passes second test\n')
                fout.write('%s\n' % (os.path.basename(fname)))

fout.close()
