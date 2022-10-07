'''
plot_cfad_npol_over_dow.py
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

#indir = '/home/disk/bob/olympex/zebra/moments/npol_qc2/rhi'
indir = '/home/disk/bob/olympex/zebra/moments/npol_qc2_ctrOnDow_to20km/rhi'
outdir = indir

#fname = 'npol_dbz_cfad_over_dow_to20km_adj-0.0.txt'
fname = 'npol_dbz_cfad_over_dow_to20km_adj-4.5.txt'
#fname = 'npol_dbz_tot_cfad_over_dow_to20km_adj-0.0db.txt'

#title = 'NPOL CZ over DOW (to 20km) Norm CFAD'
title = 'NPOL CZ-4.5dB over DOW (to 20km) Norm CFAD'
#title = 'NPOL ZZ over DOW (to 20km) Norm CFAD'

# values from input cfad text file
min_z = 0.0    # km
delta_z = 0.5  # km
min_refl = -25
max_refl = 50.0  #dBZ
delta_refl = 0.5
refl_offset = -4.5
#refl_offset = -0.0

#values for plot
min_refl_for_plot = 17.0  #dBZ
max_refl_for_plot = 45.0  #dBZ
min_z_for_plot = 1    #km
max_z_for_plot = 8    #km
min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz is 0.5km

xaxis_label = 'Reflectivity (DBZ)'
yaxis_label = 'Height (km)'

outfile = 'npol_dbz_cfad_over_dow_to20km_adj'+str(refl_offset)+'_'+str(int(min_refl_for_plot))+'db.png'
#outfile = 'npol_dbz_tot_cfad_over_dow_to20km_adj'+str(refl_offset)+'_'+str(int(min_refl_for_plot))+'db.png'

# ---------------------------END INPUTS---------------------------------

num_refl_intervals = (int)((max_refl-min_refl)/delta_refl)

#read text file
infile = indir + '/' + fname
rawdata = np.loadtxt(infile)

#use only refl >= min_refl_for_plot - input file has values from 'min_refl' to 'max_refl'
# in intervals of 'delta_refl'
start_bin =  (int)(math.floor((min_refl_for_plot - min_refl)/delta_refl))
rawdata = rawdata[:,start_bin:num_refl_intervals]

#find mode at each height starting at 2km and write to file
fid = open(indir+'/npol_modes_dbz_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt', 'w')
#fid = open(indir+'/npol_modes_dbz_tot_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt', 'w')
refl = np.zeros(rawdata.shape[0]-min_level_for_plot)
for ilevel in range(min_level_for_plot,rawdata.shape[0]):
    index = np.argmax(rawdata[ilevel,:])
    refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)
refl.tofile(fid,sep="\t",format="%f")
fid.close()

#normalize data
rawdata_norm = rawdata / np.max(rawdata)

#bin dimensions for refl
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

plt.savefig(indir+'/'+outfile, bbox_inches='tight')

