'''
plot_cfad_d3r_ka.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

indir = '/home/disk/bob/olympex/zebra/moments/d3r_ka/rhi'
outdir = indir
#title = 'D3R Ka Total Normalized CFAD (DBZ)'
#title = 'D3R Ka East Normalized CFAD (DBZ)'
#title = 'D3R Ka West Normalized CFAD (DBZ)'
title = 'D3R Ka East-West Normalized CFAD (DBZ)'
#outfile = 'd3r_ka_dbz_cfad_total.png'
#outfile = 'd3r_ka_dbz_cfad_east.png'
#outfile = 'd3r_ka_dbz_cfad_west.png'
outfile = 'd3r_ka_dbz_cfad_east-west.png'
xaxis_label = 'Reflectivity (dBZ)'
#xaxis_label = 'Differential Reflectivity (ZDR)'

#type = 'single'
type = 'diff'

if type == 'diff':
    
    #read text files
    infile1 = indir+'/d3r_ka_dbz_cfad_east.txt'
    #infile1 = indir+'/d3r_ka_zdr_cfad_east.txt'
    rawdata1 = np.loadtxt(infile1)
    infile2 = indir+'/d3r_ka_dbz_cfad_west.txt'
    #infile2 = indir+'/d3r_zdr_cfad_west.txt'
    rawdata2 = np.loadtxt(infile2)

    #normalize data
    rawdata_norm1 = rawdata1 / np.max(rawdata1)
    rawdata_norm2 = rawdata2 / np.max(rawdata2)

    #diff two input arrays
    rawdata_norm = rawdata_norm1 - rawdata_norm2

elif type == 'single':
    
    #read text file
    #infile = indir+'/d3r_ka_dbz_cfad_total.txt'
    #infile = indir+'/d3r_ka_dbz_cfad_east.txt'
    infile = indir+'/d3r_ka_dbz_cfad_west.txt'
    #infile = indir+'/d3r_ka_zdr_cfad_west.txt'
    rawdata = np.loadtxt(infile)

    #normalize data
    rawdata_norm = rawdata / np.max(rawdata)

else:

    print "Unknown type -- must be either 'diff' or 'single'"
    exit()
    
#bin dimensions
# for refl
#xcenter = np.linspace(0.5,49.5,50)
#xcenter = np.linspace(-14.5,49.5,65)
xcenter = np.linspace(-9.5,34.5,45)
ycenter = np.linspace(0.0,12,25)
# for zdr
#xcenter = np.linspace(-4.9,4.9,50)
#ycenter = np.linspace(0.0,12,25)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

if type == 'single':
    cbar_levels = np.linspace(0.05,1,20)
else:
    cbar_levels = np.linspace(-0.1,0.175,12)   #dbz
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
ax0.set_ylim(1,8)
# for refl
#ax0.set_xlim(0,50)
#ax0.set_xlim(-15,50)
ax0.set_xlim(-10,35)
# for zdr
#ax0.set_xlim(-1.5,2.5)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)
    
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'/'+outfile, bbox_inches='tight')

