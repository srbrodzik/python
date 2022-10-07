#!/usr/bin/python3

'''
plot_cfad_gpm_Dm.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot

27 Apr 2022
Modified to use netcdf file as input
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
import netCDF4 as nc4
from datetime import datetime

# use magic command to allow multiple plots
#%matplotlib

# --------------------------START INPUTS---------------------------------

# OLYMPEX cfad for
#indir = '/home/disk/bob/gpm/n1_oly_lynn_land_ku/classify/ex_data'
indir = '/home/disk/bob/gpm/nam_ku/classify/class_data_v07/stats_class_v12s/Cfad_Dm/08'
outdir = indir

apply_filter = False

#output_format = 'eps'
output_format = 'pdf'

# input file (choose one)
#fname = 'infoCfad_Dm_Conv_EchoCores_08_2018_NAM_v12s.nc'
#fname = 'infoCfad_Dm_Stra_EchoCores_08_2018_NAM_v12s.nc'
fname = 'infoCfad_Dm_ShallowIsol_EchoCores_08_2018_NAM_v12s.nc'

parts = fname.split('_')
stormType = parts[2]
month = parts[4]
year = parts[5]
dtObj = datetime.strptime(year+month,'%Y%m')
dtStr = dtObj.strftime('%b %Y')

# plot title - storm
title = 'GPMKu '+stormType+' Storm Dm - '+dtStr

# plot title - cores (uncomment one)
#coreType = 'DWC'
#coreType = 'WCC'
#coreType = 'DCC'
#coreType = 'BS'
#coreType = 'SHI'
#title = 'GPMKu '+coreType+' Core Dm - '+dtStr

# values from input cfad text file
min_z = 0.0      # km
delta_z = 0.125  # km
min_Dm = 0       # mm
max_Dm = 5.0     # mm
delta_Dm = 0.1

#values for plot
min_Dm_for_plot = 0.0 # mm
max_Dm_for_plot = 5.0 # mm
min_z_for_plot = 1    # km
max_z_for_plot = 15   # km
#min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz of smoothed data is 0.5km
min_level_for_plot = 16  # 2km - lowest height is 0km and deltaz of full data is 0.125km

xaxis_label = 'Dm (mm)'
yaxis_label = 'Height (km)'

# For Storm plots
if apply_filter:
    if output_format == 'eps':
        outfile = 'gpm-ku_cfad_Dm_'+stormType+'_Storm_'+year+month+'_filtered.eps'
    elif output_format == 'pdf':
        outfile = 'gpm-ku_cfad_Dm_'+stormType+'_Storm_'+year+month+'_filtered.pdf'
    else:
        print('Output file not eps or pdf . . . exiting.')
        sys.exit()
else:
    if output_format == 'eps':
        outfile = 'gpm-ku_cfad_Dm_'+stormType+'_Storm_'+year+month+'.eps'
    elif output_format == 'pdf':
        outfile = 'gpm-ku_cfad_Dm_'+stormType+'_Storm_'+year+month+'.pdf'
    else:
        print('Output file not eps or pdf . . . exiting.')
        sys.exit()
"""
# For Core Plots
if apply_filter:
    if output_format == 'eps':
        outfile = 'gpm-ku_cfad_Dm_'+coreType+'_Core_'+year+month+'_filtered.eps'
    elif output_format == 'pdf':
        outfile = 'gpm-ku_cfad_Dm_'+coreType+'_Core_'+year+month+'_filtered.pdf'
    else:
        print('Output file not eps or pdf . . . exiting.')
        sys.exit()
else:
    if output_format == 'eps':
        outfile = 'gpm-ku_cfad_Dm_'+coreType+'_Core_'+year+month+'.eps'
    elif output_format == 'pdf':
        outfile = 'gpm-ku_cfad_Dm_'+coreType+'_Core_'+year+month+'.pdf'
    else:
        print('Output file not eps or pdf . . . exiting.')
        sys.exit()
"""
# ---------------------------END INPUTS---------------------------------

num_Dm_intervals = (int)((max_Dm-min_Dm)/delta_Dm) + 1

# read netcdf file
ncid = nc4.Dataset(indir+'/'+fname, 'r')

# for storms
rawdata = ncid.variables['Dm_CFAD_Full'][:]

# for cores (choose one)
#rawdata = ncid.variables['Dm_CFAD_Core'][:][:]   # for Conv cores
#rawdata = ncid.variables['Dm_CFAD_Core'][:]      # for Stra or SHI cores
ncid.close()

"""
# for convective cores
if coreType == 'DWC':
    rawdata = rawdata[0,:,:]
elif coreType == 'WCC':
    rawdata = rawdata[1,:,:]
elif coreType == 'DCC':
    rawdata = rawdata[2,:,:]
else:
    print('No such coretype = ',coreType)
    sys.exit()
"""
#use only Dm >= min_Dm_for_plot - input file has values from 'min_Dm' to 'max_Dm'
# in intervals of 'delta_Dm'
start_bin =  (int)(math.floor((min_Dm_for_plot - min_Dm)/delta_Dm))
rawdata = rawdata[:,start_bin:num_Dm_intervals]
rawdata_norm = np.zeros(np.shape(rawdata))
#rawdata_smoothed = rawdata[::4,:] + rawdata[1::4,:] + rawdata[2::4,:] + rawdata[3::4,:]
#normdata_smoothed = np.zeros(np.shape(rawdata_smoothed))

#normalize data
max_count = np.max(rawdata)
rawdata_norm = rawdata / max_count

if apply_filter:
    #apply smoothing function
    sigma=1
    #rawdata_norm_smooth = gaussian(rawdata_norm,sigma)
    rawdata_norm_smooth = gaussian_filter(rawdata_norm,sigma=sigma)
        
#bin dimensions
#xcenter = np.linspace(min_Dm_for_plot+delta_Dm/2,max_Dm-delta_Dm/2,(max_Dm-min_Dm_for_plot)/delta_Dm)
xcenter = np.linspace( int(min_Dm_for_plot),
                       int(max_Dm),
                       int(((max_Dm-min_Dm_for_plot)/delta_Dm)+1) )
ycenter = np.linspace(min_z,(rawdata.shape[0]-1)*delta_z,rawdata.shape[0])

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(0.05,1,20)

fig = plt.figure()
fig.set_size_inches(4,4)
#fig.set_size_inches(2,2)
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
ax0.set_xlim(min_Dm_for_plot,max_Dm_for_plot)
ax0.text(2.5,13,'max count: '+str(max_count), fontsize=10)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)

cbar.ax.set_ylabel(r'Normalized Frequency')

if output_format == 'eps':
    plt.savefig(outdir+'/'+outfile, bbox_inches='tight', format='eps', dpi=1000)
elif output_format == 'pdf':
    plt.savefig(outdir+'/'+outfile, bbox_inches='tight', format='pdf', dpi=1000)
else:
    print('Output file not eps or pdf . . . exiting.')
    sys.exit()
