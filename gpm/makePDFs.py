import os
import gzip
import netCDF4 as nc4
import numpy as np
import shutil
import tempfile

# For plotting
%matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

# User inputs
baseDir = '/home/disk/archive3/trmm/SAM/interp_data'
year = ['1998','1999',
        '2000','2001','2002','2003','2004','2005','2006','2007','2008','2009',
        '2010','2011','2012','2013']
month = ['11','12']
outDir = '/home/disk/bob/trmm_v7/southAmer/classify/plots'
outFile_info = 'HIST_MAXREFL_CONV_INFO'
outFile_edge = 'HIST_MAXREFL_CONV_EDGE'
outFile_hist = 'HIST_MAXREFL_CONV'

CONV = 2
MIN_LAT = -36.0
MAX_LAT = -30.0
MIN_LON = -70.0
MAX_LON = -55.0

# Declare histogram
binSize = 1   # dBZ
minRefl = 0   # dBZ
maxRefl = 75  # dBZ
hist_total = np.zeros( (maxRefl-minRefl)/binSize )

for iyear in year:
    for imonth in month:
        for fname in os.listdir(baseDir+'/'+iyear+'/'+imonth):
            print fname    
            if fname.endswith(('.gz','.nc')) :
                if fname.endswith('.gz'):
                    infile = gzip.open(baseDir+'/'+iyear+'/'+imonth+'/'+fname, 'rb')
                    tmp = tempfile.NamedTemporaryFile(delete=False)
                    shutil.copyfileobj(infile, tmp)
                    infile.close()
                    tmp.close()
                    ncid = nc4.Dataset(tmp.name,'r')
                    os.unlink(tmp.name)
                else:
                    ncid = nc4.Dataset(baseDir+'/'+iyear+'/'+imonth+'/'+fname,'r')

                # Read refl, raintype, lat and lon vars
                refl = np.array(ncid.variables['corr_Zfactor'])
                (ntime,nalt,nlat,nlon) = refl.shape
                refl = np.squeeze(refl)
                refl_missing = ncid.variables['corr_Zfactor'].missing_value
            
                rt = np.array(ncid.variables['rain_type'])
                rt = np.squeeze(rt)
                rt_missing = ncid.variables['rain_type'].missing_value

                lat = np.array(ncid.variables['latitude'])
                lon = np.array(ncid.variables['longitude'])
            
                ncid.close()

                # Only consider values in lat/lon limits
                latli = np.argmin( np.abs( lat - MIN_LAT ) )
                latui = np.argmin( np.abs( lat - MAX_LAT ) )
                lonli = np.argmin( np.abs( lon - MIN_LON ) )
                lonui = np.argmin( np.abs( lon - MAX_LON ) )
                rt_sub = rt[latli:latui+1,lonli:lonui+1]
                if rt_sub.size > 0:
                    refl_sub = refl[:,latli:latui+1,lonli:lonui+1]
            
                    # Find max refl in each column
                    max_refl_sub = np.zeros(rt_sub.shape,dtype=float)+refl_missing
                    max_refl_sub = np.amax(refl_sub,axis=0)

                    # Consider only CONV columns
                    max_refl_sub_conv = max_refl_sub[ np.logical_and( max_refl_sub != refl_missing, rt_sub == CONV) ]

                    # Bin max_refl_sub_conv and add to to hist_total
                    [hist,bin_edges] = np.histogram(max_refl_sub_conv,bins=maxRefl-minRefl,range=(minRefl,maxRefl) ) 
                    hist_total = hist_total + hist

# Save hist_total and bin_edges
outInfo=np.array([minRefl,maxRefl,binSize])
np.savetxt(outDir+'/'+outFile_info,outInfo,fmt='%4.1f')
np.savetxt(outDir+'/'+outFile_edge,bin_edges,fmt='%3.1f')
np.savetxt(outDir+'/'+outFile_hist,hist_total,fmt='%d',delimiter=' ')

#plot histogram
arr = np.loadtxt(outDir+'/'+outFile_hist)
bin = np.arange(0,75)
plt.plot(bin,arr)
plt.xlabel('Reflectivity')
plt.ylabel('Counts')
plt.title('TRMM Max Refl For Conv Pixels [70W-55W,36S-30S]')
plt.grid(True)
plt.show()

np.sum(arr[0:30])
np.sum(arr[31:])
