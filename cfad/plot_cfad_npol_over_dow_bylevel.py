'''
plot_cfad_npol_over_dow_bylevel.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#indir = '/home/disk/bob/olympex/zebra/moments/npol_qc2/rhi'
indir = '/home/disk/bob/olympex/zebra/moments/npol_qc2_ctrOnDow_to20km/rhi'
outdir = indir

#fname = 'npol_dbz_cfad_over_dow.txt'
fname = 'npol_dbz_cfad_over_dow_to20km.txt'
#title = 'NPOL over DOW Normalized CFAD (DBZ) By Level'
title = 'NPOL over DOW (to 20km) Normalized CFAD (DBZ) By Level'
#outfile = 'npol_dbz_cfad_over_dow_bylevel.png'
outfile = 'npol_dbz_cfad_over_dow_to20km_bylevel.png'
xaxis_label = 'Reflectivity (DBZ)'
yaxis_label = 'Height (km)'

#read text file
infile = indir + '/' + fname
rawdata = np.loadtxt(infile)
normdata = np.zeros(np.shape(rawdata))

#normalize data
for i in range(0,rawdata.shape[0]):
    if np.max(rawdata[i,:]) != 0:
        data_level = rawdata[i,:] / np.max(rawdata[i,:]) #pdf for each level
        normdata[i,:] = data_level

#bin dimensions for refl
xcenter = np.linspace(-24.5,49.5,rawdata.shape[1])
#ycenter = np.linspace(0.25,12.25,rawdata.shape[0])
ycenter = np.linspace(0.25,9.25,rawdata.shape[0])

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(0.05,1,20)
#cbar_levels = np.linspace(0.05,0.4,20)

fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title(title)
cfad = ax0.contourf(X_grid, Y_grid, normdata, levels=cbar_levels,cmap='jet')

ax0.set_xlabel(xaxis_label)
ax0.set_ylabel(yaxis_label)
#ax0.set_ylim(1,12)
ax0.set_ylim(1,8)
ax0.set_xlim(-25,50)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=cbar_levels)
    
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(indir+'/'+outfile, bbox_inches='tight')

