'''
plot_cfad_gpm.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

outdir = '/home/disk/meso-home/jzagrod/Misc/For_Stacy/'

#read text file
infile = '/home/disk/meso-home/jzagrod/Misc/For_Stacy/cfad_total.txt'
rawdata = np.loadtxt(infile)

#normalize data
rawdata_norm = rawdata / np.max(rawdata)

#bin dimensions
xcenter = np.linspace(0.5,49.5,50)
ycenter = np.linspace(0.5,12,24)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(0.05,1,20)

fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title('Normalized CFAD')
cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel('Reflectivity (dBZ)')
ax0.set_ylabel('Height (km)')
ax0.set_ylim(0,12)
ax0.set_xlim(0,50)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=np.linspace(0.05,1,20))
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'cfad_test.png', bbox_inches='tight')

