'''
plot_cfad_gpm_diff_lynn_bylevel_ka.py
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

# use magic command to allow multiple plots
#%matplotlib

# --------------------------START INPUTS---------------------------------

indir_land = '/home/disk/bob/gpm/n1_oly_lynn_land_ka_hs/classify/ex_data'
indir_ocean = '/home/disk/bob/gpm/n1_oly_lynn_ocean_ka_hs/classify/ex_data'
outdir = indir_land
#outdir = '.'

apply_filter = True

fname_land = 'gpm-ka_cfad_lynn_land_total_NEW.txt'
fname_ocean = 'gpm-ka_cfad_lynn_ocean_total_NEW.txt'

title = 'GPMKaHS Land-Ocean Norm By Level CFAD'

# values from input cfad text file
min_z = 0.0    # km
delta_z = 0.25  # km
min_refl = -25
max_refl = 50.0  #dBZ
delta_refl = 0.5

#values for plot
min_refl_for_plot = 10.0  #dBZ
max_refl_for_plot = 45.0  #dBZ
min_z_for_plot = 1    #km
max_z_for_plot = 8    #km
#min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz of smoothed data is 0.5km
min_level_for_plot = 8   # 2km - lowest height is 0km and deltaz of full data is 0.25km

xaxis_label = 'Reflectivity (DBZ)'
yaxis_label = 'Height (km)'

if apply_filter:
    outfile = 'gpm-ka_cfad_lynn_diff_total_'+str(int(min_refl_for_plot))+'db_bylevel_filtered.eps'
else:
    outfile = 'gpm-ka_cfad_lynn_diff_total_'+str(int(min_refl_for_plot))+'db_bylevel.eps'

# ---------------------------END INPUTS---------------------------------

num_refl_intervals = (int)((max_refl-min_refl)/delta_refl)

#read text files
infile_land = indir_land+'/'+fname_land
rawdata_land = np.loadtxt(infile_land)
infile_ocean = indir_ocean+'/'+fname_ocean
rawdata_ocean = np.loadtxt(infile_ocean)

#use only refl >= min_refl_for_plot - input file has values from 'min_refl' to 'max_refl'
# in intervals of 'delta_refl'
start_bin =  (int)(math.floor((min_refl_for_plot - min_refl)/delta_refl))
rawdata_land = rawdata_land[:,start_bin:num_refl_intervals]
rawdata_ocean = rawdata_ocean[:,start_bin:num_refl_intervals]
#rawdata_smoothed = rawdata[::4,:] + rawdata[1::4,:] + rawdata[2::4,:] + rawdata[3::4,:]
#normdata_smoothed = np.zeros(np.shape(rawdata_smoothed))

#output smoothed data for testing only
#numSmoothedLevels = rawdata_smoothed.shape[0]
#fid_smoothed = open(outdir+'/gpm-ka-over-land_cfad_total_smoothed_'+str(int(min_refl_for_plot))+'db.txt','w')
#for ilevel in range(0,numSmoothedLevels):
#    rawdata_smoothed[ilevel,:].tofile(fid_smoothed,sep="\t",format="%d")
#    fid_smoothed.write('\n')
#fid_smoothed.close()

#find mode at each height starting at min_level_for_plot (2km) and write to file
#fid = open(indir+'/gpm_modes_over_land_adj_all_levels_'+str(int(min_refl_for_plot))+'db.txt', 'w')
#refl = np.zeros(rawdata.shape[0]-min_level_for_plot)
#for ilevel in range(min_level_for_plot,rawdata.shape[0]):
#    index = np.argmax(rawdata[ilevel,:])
#    refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)

#fid = open(indir+'/gpm_modes_over_land_adj_smoothed_'+str(int(min_refl_for_plot))+'db.txt', 'w')
#fid = open(indir+'/gpm_modes_over_land_raw_smoothed_'+str(int(min_refl_for_plot))+'db.txt', 'w')
#refl = np.zeros(rawdata_smoothed.shape[0]-min_level_for_plot)
#for ilevel in range(min_level_for_plot,rawdata_smoothed.shape[0]):
#    index = np.argmax(rawdata_smoothed[ilevel,:])
#    refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)

#refl.tofile(fid,sep="\t",format="%f")
#fid.close()

#normalize data
#rawdata_norm = rawdata / np.max(rawdata)

#normalize data by level
rawdata_land_norm = np.zeros(np.shape(rawdata_land))
for i in range(0,rawdata_land.shape[0]):
    if np.max(rawdata_land[i,:]) != 0:
        data_level = rawdata_land[i,:] / np.max(rawdata_land[i,:]) #pdf for each level
        rawdata_land_norm[i,:] = data_level
rawdata_ocean_norm = np.zeros(np.shape(rawdata_ocean))
for i in range(0,rawdata_ocean.shape[0]):
    if np.max(rawdata_ocean[i,:]) != 0:
        data_level = rawdata_ocean[i,:] / np.max(rawdata_ocean[i,:]) #pdf for each level
        rawdata_ocean_norm[i,:] = data_level

#diff the two input arrays
rawdata_diff_norm = rawdata_land_norm - rawdata_ocean_norm

if apply_filter:
    #apply smoothing function
    sigma=1
    #rawdata_diff_norm_smooth = gaussian(rawdata_diff_norm,sigma)
    rawdata_diff_norm_smooth = gaussian_filter(rawdata_diff_norm,sigma=sigma)

#replace all zeros with NaN's to pretty up plot
if apply_filter:
    rawdata_diff_norm_smooth[rawdata_diff_norm_smooth==0]=np.NaN
else:
    rawdata_diff_norm[rawdata_diff_norm==0]=np.NaN
    
#bin dimensions
xcenter = np.linspace(min_refl_for_plot+delta_refl/2,max_refl-delta_refl/2,(max_refl-min_refl_for_plot)/delta_refl)
ycenter = np.linspace(min_z,(rawdata_land.shape[0]-1)*delta_z,rawdata_land.shape[0])

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(-0.95,0.95,20)
#cbar_levels = np.linspace(-1.00,1.00,20)   #dbz above 2km
#cbar_levels = np.linspace(-0.35,0.35,14)   #dbz above 2km

fig = plt.figure()
fig.set_size_inches(4,4)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title(title)
if apply_filter:
    #cfad = ax0.contourf(X_grid, Y_grid, rawdata_diff_norm_smooth, levels=cbar_levels,cmap='jet')
    cfad = ax0.contourf(X_grid, Y_grid, rawdata_diff_norm_smooth, levels=cbar_levels,cmap='seismic')
else:
    cfad = ax0.contourf(X_grid, Y_grid, rawdata_diff_norm, levels=cbar_levels,cmap='seismic')
    
ax0.set_xlabel(xaxis_label)
ax0.set_ylabel(yaxis_label)
ax0.set_ylim(min_z_for_plot,max_z_for_plot)
ax0.set_xlim(min_refl_for_plot,max_refl_for_plot)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)

cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'/'+outfile, bbox_inches='tight', format='eps', dpi=1000)

plt.savefig(outdir+'/'+outfile, bbox_inches='tight', format='jpg', dpi=1000)

