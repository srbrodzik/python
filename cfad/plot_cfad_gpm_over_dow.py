'''
plot_cfad_gpm.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot
'''

import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

indir = '/home/disk/bob/gpm/n1_oly_npol_ku/classify/ex_data'
#indir = '/home/disk/bob/gpm/n1_oly_npol_ku/classify/ex_data_west'
outdir = indir

fname = 'gpm-ku-over-dow_cfad_total.txt'
#fname = 'gpm-ku-over-npol-west_cfad_total.txt'

title = 'GPM-Ku Total Normalized CFAD - Over DOW'
#title = 'GPM-Ku Total Normalized CFAD - Over NPOL - West'

outfile = 'gpm-ku-over-dow_cfad_total.png'

min_refl = -25
max_refl = 50
interval = 0.5
numIntervals = (int)((max_refl-min_refl)/interval)

min_refl_for_plot = 10  # dBZ
min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz of smoothed data is 0.5km

xaxis_label = 'Reflectivity (DBZ)'
yaxis_label = 'Height (km)'

#read text file
infile = indir+'/'+fname
rawdata = np.loadtxt(infile)

#use only refl >= min_dbz_for_plot - input file has values from 'min_refl' to 'max_refl'
# in intervals of 'interval'
start_bin =  (int)(math.floor((min_refl_for_plot - min_refl)/interval))
rawdata = rawdata[:,start_bin:numIntervals]
rawdata_smoothed = rawdata[::4,:] + rawdata[1::4,:] + rawdata[2::4,:] + rawdata[3::4,:]

#find mode at each height starting at 2km and write to file
fid = open(indir+'/gpm_modes.txt', 'w')
refl = np.zeros(rawdata_smoothed.shape[0]-min_level_for_plot)
for ilevel in range(min_level_for_plot,rawdata_smoothed.shape[0]):
    index = np.argmax(rawdata_smoothed[ilevel,:])
    refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(interval/2.)) + (index*interval)
refl.tofile(fid,sep="\t",format="%f")
fid.close()

#normalize data
rawdata_norm = rawdata / np.max(rawdata)

#bin dimensions
xcenter = np.linspace(9.5,49.5,80)  # for 10 to 50 refl range in 0.5 increments
#xcenter = np.linspace(0.5,49.5,50)
#xcenter = np.linspace(-24.5,49.5,75)
#ycenter = np.linspace(0.25,9.25,19)
ycenter = np.linspace(0.0,21.875,176)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(0.05,1,20)

fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title(title)
cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel(xaxis_label)
ax0.set_ylabel(yaxis_label)
ax0.set_ylim(1,8)
#ax0.set_xlim(-25,50)
ax0.set_xlim(10,40)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)

cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'/'+outfile, bbox_inches='tight')

