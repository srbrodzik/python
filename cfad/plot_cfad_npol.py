'''
plot_cfad_npol.py
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

# values from input cfad text file
min_z = 0.0    # km
delta_z = 0.5  # km
min_refl = -25
max_refl = 50.0  #dBZ
delta_refl = 0.5
refl_offset = -4.5
#refl_offset = -0.0

#values for plot
min_refl_for_plot = -10.0  #dBZ
max_refl_for_plot = 40.0  #dBZ
min_z_for_plot = 1    #km
max_z_for_plot = 8    #km
min_level_for_plot = 4  # 2km - lowest height is 0km and deltaz is 0.5km
xaxis_label = 'Reflectivity (dBZ)'
#xaxis_label = 'Differential Reflectivity (ZDR)'
yaxis_label = 'Height (km)'
#title = 'NPOL CZ'+str(refl_offset)+'dB Total Norm CFAD'
#title = 'NPOL CZ'+str(refl_offset)+'dB East Norm CFAD'
#title = 'NPOL CZ'+str(refl_offset)+'dB West Norm CFAD'
title = 'NPOL CZ'+str(refl_offset)+'dB East-West Norm CFAD'

indir = '/home/disk/bob/olympex/zebra/moments/npol_qc2/rhi'
outdir = indir

fname_total = 'npol_dbz_cfad_total_adj'+str(refl_offset)+'db.txt'
fname_east = 'npol_dbz_cfad_east_adj'+str(refl_offset)+'db.txt'
fname_west = 'npol_dbz_cfad_west_adj'+str(refl_offset)+'db.txt'

#outfile = 'npol_dbz_cfad_total_adj'+str(refl_offset)+'db.png'
#outfile = 'npol_dbz_cfad_east_adj'+str(refl_offset)+'db.png'
#outfile = 'npol_dbz_cfad_west_adj'+str(refl_offset)+'db.png'
outfile = 'npol_dbz_cfad_east-west_adj'+str(refl_offset)+'db.png'

#type = 'single'
type = 'diff'

# ---------------------------END INPUTS---------------------------------

num_refl_intervals = (int)((max_refl-min_refl)/delta_refl)

#read text files
if type == 'diff':
    infile1 = indir+'/'+fname_east
    rawdata1 = np.loadtxt(infile1)
    infile2 = indir+'/'+fname_west
    rawdata2 = np.loadtxt(infile2)
elif type == 'single':
    #infile = indir+'/'+fname_total
    #infile = indir+'/'+fname_east
    infile = indir+'/'+fname_west
    rawdata = np.loadtxt(infile)
else:
    print "Unknown type -- must be either 'diff' or 'single'"
    exit()
    
#use only refl >= min_refl_for_plot - input file(s) has values from 'min_refl' to 'max_refl'
# in intervals of 'delta_refl'
start_bin =  (int)(math.floor((min_refl_for_plot - min_refl)/delta_refl))
if type == 'diff':
    rawdata1 = rawdata1[:,start_bin:num_refl_intervals]
    rawdata2 = rawdata2[:,start_bin:num_refl_intervals]
elif type == 'single':
    rawdata = rawdata[:,start_bin:num_refl_intervals]
  
#find mode at each height starting at min_level_for_plot (2km) and write to file
if type == 'diff':
    fid1 = open(indir+'/npol_modes_east_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt', 'w')
    refl = np.zeros(rawdata1.shape[0]-min_level_for_plot)
    for ilevel in range(min_level_for_plot,rawdata1.shape[0]):
        index = np.argmax(rawdata1[ilevel,:])
        refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)
    refl.tofile(fid1,sep="\t",format="%f")
    fid1.close()
    fid2 = open(indir+'/npol_modes_west_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt', 'w')
    refl = np.zeros(rawdata2.shape[0]-min_level_for_plot)
    for ilevel in range(min_level_for_plot,rawdata2.shape[0]):
        index = np.argmax(rawdata2[ilevel,:])
        refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)
    refl.tofile(fid2,sep="\t",format="%f")
    fid2.close()
else:
    #fid = open(indir+'/npol_modes_total_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt', 'w')
    #fid = open(indir+'/npol_modes_east_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt', 'w')
    fid = open(indir+'/npol_modes_west_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt', 'w')
    refl = np.zeros(rawdata.shape[0]-min_level_for_plot)
    for ilevel in range(min_level_for_plot,rawdata.shape[0]):
        index = np.argmax(rawdata[ilevel,:])
        refl[ilevel-min_level_for_plot] = (min_refl_for_plot+(delta_refl/2.)) + (index*delta_refl)
    refl.tofile(fid,sep="\t",format="%f")
    fid.close()

#normalize data
if type == 'diff':
    rawdata_norm1 = rawdata1 / np.max(rawdata1)
    rawdata_norm2 = rawdata2 / np.max(rawdata2)
    rawdata_norm = rawdata_norm1 - rawdata_norm2
else:
    rawdata_norm = rawdata / np.max(rawdata)
    
# bin dimensions for refl
xcenter = np.linspace(min_refl_for_plot+delta_refl/2,max_refl-delta_refl/2,(max_refl-min_refl_for_plot)/delta_refl)
if type == 'diff':
    ycenter = np.linspace(min_z,(rawdata1.shape[0]-1)*delta_z,rawdata1.shape[0])
else:
    ycenter = np.linspace(min_z,(rawdata.shape[0]-1)*delta_z,rawdata.shape[0])
## for zdr
#xcenter = np.linspace(-4.9,4.9,50)
#ycenter = np.linspace(0.5,12,24)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

if type == 'diff':
    cbar_levels = np.linspace(-0.3,0.3,14)   #dbz above 2km
    #cbar_levels = np.linspace(-0.6,0.3,14)   #dbz to ground
    #cbar_levels = np.linspace(-1.0,0.1,12)   #zdr
else:
    cbar_levels = np.linspace(0.05,1,20)
    #cbar_levels_to_print = np.linspace(0.05,0.95,10)

fig = plt.figure()
#fig.set_size_inches(6,6)
fig.set_size_inches(4,4)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title(title)
cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel(xaxis_label)
ax0.set_ylabel(yaxis_label)
ax0.set_ylim(min_z_for_plot,max_z_for_plot)
ax0.set_xlim(min_refl_for_plot,max_refl_for_plot)  #refl
#ax0.set_xlim(min_zdr_for_plot,max_zdr_for_plot)   #zdr
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)
#cbar = plt.colorbar(cfad,cax,ticks=cbar_levels_to_print)
   
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'/'+outfile, bbox_inches='tight')

