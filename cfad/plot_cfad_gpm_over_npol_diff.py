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

#outdir = '/home/disk/bob/gpm/n1_oly_npol_ku/classify/ex_data/'
outdir = '/home/disk/bob/gpm/n1_oly_npol_ku/classify/ex_data/'

#read text file
infile1 = '/home/disk/bob/gpm/n1_oly_npol_ku/classify/ex_data_east/gpm-ku-over-npol-east_cfad_total.txt'
rawdata1 = np.loadtxt(infile1)
infile2 = '/home/disk/bob/gpm/n1_oly_npol_ku/classify/ex_data_west/gpm-ku-over-npol-west_cfad_total.txt'
rawdata2 = np.loadtxt(infile2)

#normalize data
rawdata_norm1 = rawdata1 / np.max(rawdata1)
rawdata_norm2 = rawdata2 / np.max(rawdata2)

#diff two input arrays
rawdata_norm = rawdata_norm1 - rawdata_norm2

#bin dimensions
xcenter = np.linspace(0.5,49.5,50)
ycenter = np.linspace(0.0,21.875,176)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(-0.6,0.3,14)   #dbz

fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title('GPM-Ku Total Normalized CFAD - Over NPOL - East-West')
cfad = ax0.contourf(X_grid, Y_grid, rawdata_norm, levels=cbar_levels,cmap='jet')

ax0.set_xlabel('Reflectivity (dBZ)')
ax0.set_ylabel('Height (km)')
ax0.set_ylim(0,8)
ax0.set_xlim(0,50)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=np.linspace(-0.6,0.3,14))  #dbz
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'gpm-ku-over-npol_east-west_cfad_total.png', bbox_inches='tight')

