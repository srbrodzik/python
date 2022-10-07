'''
plot_cfad_dow.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

indir = '/home/disk/bob/olympex/zebra/moments/dow_lo_qc/rhi/'
#indir = '/home/disk/bob/olympex/zebra/moments/dow_lo_qc_to20km/rhi/'
fname = 'dow_dbzhc_f_cfad_for_comp_with_npol.txt'
#fname = 'dow_dbzhc_f_cfad_for_comp_with_npol_to20km.txt'
title = 'DOW Normalized CFAD (DBZHC_F)'
#title = 'DOW (to 20km) Normalized CFAD (DBZHC_F)'
outfile = 'dow_dbzhc_f_cfad_for_comp_with_npol.png'
#outfile = 'dow_dbzhc_f_cfad_for_comp_with_npol_to20km.png'
xaxis_label = 'Reflectivity (DBZHC_F)'
yaxis_label = 'Height (km)'

#read text file
infile = indir + '/' + fname
rawdata = np.loadtxt(infile)

#normalize data
rawdata_norm = rawdata / np.max(rawdata)
#rawdata_norm = rawdata

#bin dimensions for refl
xcenter = np.linspace(-24.5,49.5,75)
ycenter = np.linspace(0.25,12.25,25)
#ycenter = np.linspace(0.25,9.25,19)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(0.05,1,20)
#cbar_levels = np.linspace(0,200000,20)

fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title(title)
cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel(xaxis_label)
ax0.set_ylabel(yaxis_label)
#ax0.set_ylim(1,12)
ax0.set_ylim(1,8)
ax0.set_xlim(-25,50)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)
    
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(indir+outfile, bbox_inches='tight')

