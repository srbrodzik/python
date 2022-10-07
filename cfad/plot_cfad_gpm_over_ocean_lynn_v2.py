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
#from scipy.signal import gaussian
from scipy.ndimage.filters import gaussian_filter
import sys

# use magic command to allow multiple plots
#%matplotlib

# --------------------------START INPUTS---------------------------------

# OLYMPEX cfad for paper
indir = '/home/disk/bob/gpm/n1_oly_lynn_ocean_ku/classify/ex_data'

#indir = '/home/disk/bob/gpm/vanc_isl_lynn_ocean_ku/classify/ex_data'
#indir = '/home/disk/bob/gpm/north_ore_lynn_ocean_ku/classify/ex_data'
#indir = '/home/disk/bob/gpm/south_ore_lynn_ocean_ku/classify/ex_data'
#indir = '/home/disk/bob/gpm/north_ca_lynn_ocean_ku/classify/ex_data'
#indir = '/home/disk/bob/gpm/cent_ca_lynn_ocean_ku/classify/ex_data'
#indir = '/home/disk/bob/gpm/south_ca_lynn_ocean_ku/classify/ex_data'
outdir = indir

apply_filter = True

#output_format = 'eps'
output_format = 'pdf'

# OLYMPEX cfad for paper
#fname = 'gpm-ku_cfad_lynn_ocean_total.txt'
fname = 'gpm-ku_cfad_lynn_ocean_total_NEW_2014-18.txt'

#fname = 'gpm-ku_cfad_lynn_ocean_total_NEW.txt'

# OLYMPEX cfad for paper
title = 'GPMKu (Oly) Over Ocean Norm CFAD'

#title = 'GPMKu (VancIsl) Over Ocean Norm CFAD'
#title = 'GPMKu (NOre) Over Ocean Norm CFAD'
#title = 'GPMKu (SOre) Over Ocean Norm CFAD'
#title = 'GPMKu (NCa) Over Ocean Norm CFAD'
#title = 'GPMKu (CCa) Over Ocean Norm CFAD'
#title = 'GPMKu (SCa) Over Ocean Norm CFAD'

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
#min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz of smoothed data is 0.5km
min_level_for_plot = 16  # 2km - lowest height is 0km and deltaz of full data is 0.125km

xaxis_label = 'Reflectivity (DBZ)'
yaxis_label = 'Height (km)'

if apply_filter:
    if output_format == 'eps':
        outfile = 'gpm-ku_cfad_lynn_ocean_total_'+str(int(min_refl_for_plot))+'db_filtered.eps'
    elif output_format == 'pdf':
        outfile = 'gpm-ku_cfad_lynn_ocean_total_'+str(int(min_refl_for_plot))+'db_filtered.pdf'
    else:
        print 'Output file not eps or pdf . . . exiting.'
        sys.exit()
else:
    if output_format == 'eps':
        outfile = 'gpm-ku_cfad_lynn_ocean_total_'+str(int(min_refl_for_plot))+'db.eps'
    elif output_format == 'pdf':
        outfile = 'gpm-ku_cfad_lynn_ocean_total_'+str(int(min_refl_for_plot))+'db.pdf'
    else:
        print 'Output file not eps or pdf . . . exiting.'
        sys.exit()

# ---------------------------END INPUTS---------------------------------

num_refl_intervals = (int)((max_refl-min_refl)/delta_refl)

#read text file
infile = indir+'/'+fname
rawdata = np.loadtxt(infile)

#use only refl >= min_dbz_for_plot - input file has values from 'min_refl' to 'max_refl'
# in intervals of 'delta_refl'
start_bin =  (int)(math.floor((min_refl_for_plot - min_refl)/delta_refl))
rawdata = rawdata[:,start_bin:num_refl_intervals]
rawdata_norm = np.zeros(np.shape(rawdata))
#rawdata_smoothed = rawdata[::4,:] + rawdata[1::4,:] + rawdata[2::4,:] + rawdata[3::4,:]
#normdata_smoothed = np.zeros(np.shape(rawdata_smoothed))

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
#fid = open(indir+'/gpm_modes_over_ocean_raw_smoothed_'+str(int(min_refl_for_plot))+'db.txt', 'w')
#refl = np.zeros(rawdata_smoothed.shape[0]-min_level_for_plot)
#for ilevel in range(min_level_for_plot,rawdata_smoothed.shape[0]):
#    index = np.argmax(rawdata_smoothed[ilevel,:])
#    refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)

#refl.tofile(fid,sep="\t",format="%f")
#fid.close()

#normalize data
rawdata_norm = rawdata / np.max(rawdata)

if apply_filter:
    #apply smoothing function
    sigma=1
    #rawdata_norm_smooth = gaussian(rawdata_norm,sigma)
    rawdata_norm_smooth = gaussian_filter(rawdata_norm,sigma=sigma)
        
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
if apply_filter:
    cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm_smooth, levels=cbar_levels,cmap='jet')
else:
    cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel(xaxis_label)
ax0.set_ylabel(yaxis_label)
ax0.set_ylim(min_z_for_plot,max_z_for_plot)
ax0.set_xlim(min_refl_for_plot,max_refl_for_plot)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)

cbar.ax.set_ylabel(r'Frequency')

if output_format == 'eps':
    plt.savefig(outdir+'/'+outfile, bbox_inches='tight', format='eps', dpi=1000)
elif output_format == 'pdf':
    plt.savefig(outdir+'/'+outfile, bbox_inches='tight', format='pdf', dpi=1000)
else:
    print 'Output file not eps or pdf . . . exiting.'
    sys.exit()

