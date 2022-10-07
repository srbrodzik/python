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
import pdb

outdir = '/home/disk/meso-home/jzagrod/Misc/For_Stacy/'

#read text file
infile = '/home/disk/meso-home/jzagrod/Misc/For_Stacy/cfad_total.txt'
rawdata = np.loadtxt(infile)
normdata = np.zeros(np.shape(rawdata))

#normalize data by level
for i in range(0,24):
    data_level = rawdata[i,:] / np.sum(rawdata[i,:]) #pdf for each level
    normdata[i,:] = data_level
   

#bin dimensions
xcenter = np.linspace(0.5,49.5,50)
ycenter = np.linspace(0.5,12,24)

X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

cbar_levels = np.linspace(0.005,0.2,14)

fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title('Normalized CFAD (by level)')
cfad = ax0.contourf(X_grid, Y_grid, normdata, levels=cbar_levels,cmap='jet')

ax0.set_xlabel('Reflectivity (dBZ)')
ax0.set_ylabel('Height (km)')
ax0.set_ylim(0,12)
ax0.set_xlim(0,50)
ax0.grid()

cbar = plt.colorbar(cfad,cax,ticks=np.linspace(0.005,0.2,14))
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'cfad_test_bylevel.png', bbox_inches='tight')

pdb.set_trace()

