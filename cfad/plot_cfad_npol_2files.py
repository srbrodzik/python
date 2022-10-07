'''
plot_cfad_npol.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

outdir = '/home/disk/bob/olympex/zebra/moments/npol_to_ground/rhi/'
#title = 'NPOL Total Normalized CFAD (DBZ)'
title = 'NPOL East-West Normalized CFAD (DBZ)'
#outfile = 'npol_dbz_cfad_total.png'
outfile = 'npol_dbz_cfad_east-west.png'
xaxis_label = 'Reflectivity (dBZ)'
#xaxis_label = 'Differential Reflectivity (ZDR)'

#type = 'single'
type = 'diff'

if type == 'diff':
    
    #read text files
    infile1e = '/home/disk/bob/olympex/zebra/moments/npol_to_ground/rhi/npol1_dbz_cfad_east.txt'
    rawdata1e = np.loadtxt(infile1e)
    infile2e = '/home/disk/bob/olympex/zebra/moments/npol_to_ground/rhi/npol2_dbz_cfad_east.txt'
    rawdata2e = np.loadtxt(infile2e)
    rawdata1 = rawdata1e+rawdata2e

    infile1w = '/home/disk/bob/olympex/zebra/moments/npol_to_ground/rhi/npol1_dbz_cfad_west.txt'
    rawdata1w = np.loadtxt(infile1w)
    infile2w = '/home/disk/bob/olympex/zebra/moments/npol_to_ground/rhi/npol2_dbz_cfad_west.txt'
    rawdata2w = np.loadtxt(infile2w)
    rawdata2 = rawdata1w+rawdata2w

    #normalize data
    rawdata_norm1 = rawdata1 / np.max(rawdata1)
    rawdata_norm2 = rawdata2 / np.max(rawdata2)

    #diff two input arrays
    rawdata_norm = rawdata_norm1 - rawdata_norm2

elif type == 'single':
    
    #read text file
    infile1 = '/home/disk/bob/olympex/zebra/moments/npol_to_ground/rhi/npol1_dbz_cfad_total.txt'
    infile2 = '/home/disk/bob/olympex/zebra/moments/npol_to_ground/rhi/npol2_dbz_cfad_total.txt'
    rawdata1 = np.loadtxt(infile1)
    rawdata2 = np.loadtxt(infile2)
    rawdata = rawdata1+rawdata2

    #zero out part of array
    #rawdata[:,0:10] = 0
    
    #normalize data
    rawdata_norm = rawdata / np.max(rawdata)

else:

    print "Unknown type -- must be either 'diff' or 'single'"
    exit()
    
#bin dimensions
# for refl
xcenter = np.linspace(0.5,49.5,50)
#ycenter = np.linspace(0.5,12,24)
ycenter = np.linspace(0.0,12,25)
# for zdr
#xcenter = np.linspace(-4.9,4.9,50)
#ycenter = np.linspace(0.5,12,24)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

if type == 'single':
    cbar_levels = np.linspace(0.05,1,20)
else:
    cbar_levels = np.linspace(-0.6,0.3,14)   #dbz
    #cbar_levels = np.linspace(-1.0,0.1,12)    #zdr

fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title(title)
cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel(xaxis_label)
ax0.set_ylabel('Height (km)')
ax0.set_ylim(0,8)
# for refl
ax0.set_xlim(0,50)
# for zdr
#ax0.set_xlim(-1.5,2.5)
ax0.grid()

if type == 'single':
    cbar = plt.colorbar(cfad,cax,ticks=np.linspace(0.05,1,20))
else:
    cbar = plt.colorbar(cfad,cax,ticks=np.linspace(-0.6,0.3,14))  #dbz
    #cbar = plt.colorbar(cfad,cax,ticks=np.linspace(-1.0,0.1,12))   #zdr
    
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+outfile, bbox_inches='tight')

