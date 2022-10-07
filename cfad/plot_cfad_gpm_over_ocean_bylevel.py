'''
plot_cfad_gpm_over_ocean.py
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

# --------------------------START INPUTS---------------------------------

indir = '/home/disk/bob/gpm/n1_oly_ocean_ku/classify/ex_data'
outdir = indir

#fname = 'gpm-ku-over-ocean_cfad_total_adj.txt'
fname = 'gpm-ku-over-ocean_cfad_total.txt'

#title = 'GPMKu+1.7dB Over Ocean Norm CFAD'
title = 'GPMKu Over Ocean Norm By Level CFAD'

# values from input cfad text file
min_z = 0.0    # km
delta_z = 0.125  # km
min_refl = -25
max_refl = 50.0  #dBZ
delta_refl = 0.5

#values for plot
min_refl_for_plot = 10.0  #dBZ
max_refl_for_plot = 45.0  #dBZ
min_z_for_plot = 1    #km
max_z_for_plot = 8    #km
min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz of smoothed data is 0.5km
#min_level_for_plot = 16  # 2km - lowest height is 0km and deltaz of full data is 0.125km

xaxis_label = 'Reflectivity (DBZ)'
yaxis_label = 'Height (km)'

#outfile = 'gpm-ku-over-ocean_cfad_total_adj_'+str(int(min_refl_for_plot))+'db_bylevel.png'
outfile = 'gpm-ku-over-ocean_cfad_total_raw_'+str(int(min_refl_for_plot))+'db_bylevel.png'

# ---------------------------END INPUTS---------------------------------

num_refl_intervals = (int)((max_refl-min_refl)/delta_refl)

#read text file
infile = indir+'/'+fname
rawdata = np.loadtxt(infile)

#use only refl >= min_dbz_for_plot - input file has values from 'min_refl' to 'max_refl'
# in intervals of 'delta_refl'
start_bin =  (int)(math.floor((min_refl_for_plot - min_refl)/delta_refl))
rawdata = rawdata[:,start_bin:num_refl_intervals]
normdata = np.zeros(np.shape(rawdata))
rawdata_smoothed = rawdata[::4,:] + rawdata[1::4,:] + rawdata[2::4,:] + rawdata[3::4,:]
normdata_smoothed = np.zeros(np.shape(rawdata_smoothed))

#output smoothed data for testing only
#numSmoothedLevels = rawdata_smoothed.shape[0]
#fid_smoothed = open(outdir+'/gpm-ku-over-ocean_cfad_total_smoothed_'+str(int(min_refl_for_plot))+'db.txt','w')
#for ilevel in range(0,numSmoothedLevels):
#    rawdata_smoothed[ilevel,:].tofile(fid_smoothed,sep="\t",format="%d")
#    fid_smoothed.write('\n')
#fid_smoothed.close()

#find mode at each height starting at min_level_for_plot (2km) and write to file
#fid = open(indir+'/gpm_modes_over_ocean_adj_all_levels_'+str(int(min_refl_for_plot))+'db.txt', 'w')
#refl = np.zeros(rawdata.shape[0]-min_level_for_plot)
#for ilevel in range(min_level_for_plot,rawdata.shape[0]):
#    index = np.argmax(rawdata[ilevel,:])
#    refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)

#fid = open(indir+'/gpm_modes_over_ocean_adj_smoothed_'+str(int(min_refl_for_plot))+'db.txt', 'w')
fid = open(indir+'/gpm_modes_over_ocean_raw_smoothed_'+str(int(min_refl_for_plot))+'db.txt', 'w')
refl = np.zeros(rawdata_smoothed.shape[0]-min_level_for_plot)
for ilevel in range(min_level_for_plot,rawdata_smoothed.shape[0]):
    index = np.argmax(rawdata_smoothed[ilevel,:])
    refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)

refl.tofile(fid,sep="\t",format="%f")
fid.close()

#normalize data
#rawdata_norm = rawdata / np.max(rawdata)

#normalize data by level
for i in range(0,rawdata.shape[0]):
    if np.max(rawdata[i,:]) != 0:
        data_level = rawdata[i,:] / np.max(rawdata[i,:]) #pdf for each level
        normdata[i,:] = data_level

#bin dimensions
xcenter = np.linspace(min_refl_for_plot+delta_refl/2,max_refl-delta_refl/2,(max_refl-min_refl_for_plot)/delta_refl)
ycenter = np.linspace(min_z,(rawdata.shape[0]-1)*delta_z,rawdata.shape[0])

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(0.05,1,20)

fig = plt.figure()
fig.set_size_inches(4,4)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title(title)
cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel(xaxis_label)
ax0.set_ylabel(yaxis_label)
ax0.set_ylim(min_z_for_plot,max_z_for_plot)
ax0.set_xlim(min_refl_for_plot,max_refl_for_plot)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)

cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'/'+outfile, bbox_inches='tight')

