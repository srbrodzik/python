'''
plot_cfad_dow.py
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

#indir = '/home/disk/bob/olympex/zebra/moments/dow_lo_qc/rhi'
indir = '/home/disk/bob/olympex/zebra/moments/dow_lo_qc_to20km/rhi'
outdir = indir

#fname = 'dow_dbzhc_f_cfad_for_comp_with_npol.txt'
fname = 'dow_dbzhc_f_cfad_for_comp_with_npol_to20km.txt'

#title = 'DOW DBZHC_F Norm CFAD'
title = 'DOW DBZHC_F (to 20km) Norm CFAD'

# values from input cfad text file
min_z = 0.0    # km
delta_z = 0.5  # km
min_refl = -25
max_refl = 50.0  #dBZ
delta_refl = 0.5

#values for plot
min_refl_for_plot = 17.0  #dBZ
max_refl_for_plot = 45.0  #dBZ
min_z_for_plot = 1    #km
max_z_for_plot = 8    #km
min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz is 0.5km

xaxis_label = 'Reflectivity (DBZHC_F)'
yaxis_label = 'Height (km)'

#outfile = 'dow_dbzhc_f_cfad_for_comp_with_npol_'+str(int(min_refl_for_plot))+'db.png'
outfile = 'dow_dbzhc_f_cfad_for_comp_with_npol_to20km_'+str(int(min_refl_for_plot))+'db.png'

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
fid = open(indir+'/dow_modes_dbzhc_f_'+str(int(min_refl_for_plot))+'db.txt', 'w')
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

#xcenter = np.linspace(-24.5,49.5,150)  # for -25 to 50 refl range
#xcenter = np.linspace(10.25,49.75,80)  # for 10 to 50 refl range in 0.5 increments
#xcenter = np.linspace(15.25,49.75,70)  # for 15 to 50 refl range in 0.5 increments
#ycenter = np.linspace(0.0,12.0,25)
#ycenter = np.linspace(0.0,9.0,19)

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
#ax0.set_ylim(1,12)
#ax0.set_ylim(1,8)
#ax0.set_xlim(-25,50)  # for -25 to 50 refl range
#ax0.set_xlim(10,40)   # for 10 to 40 refl range
#ax0.set_xlim(15,40)   # for 15 to 40 refl range
ax0.set_ylim(min_z_for_plot,max_z_for_plot)
ax0.set_xlim(min_refl_for_plot,max_refl_for_plot)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)
    
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'/'+outfile, bbox_inches='tight')

