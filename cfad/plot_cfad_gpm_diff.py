'''
plot_cfad_gpm_diff.py
Joe Zagrodnik 28-Sept-2016
Read Stacy's text file
Make CFAD plot
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

outdir = '/home/disk/bob/gpm/n1_oly_ocean_ku/classify/ex_data/'

#read text file
inFileKa = '/home/disk/bob/gpm/n1_oly_ocean_ka_hs/classify/ex_data/gpm-ka_cfad_total.txt'
inFileKu = '/home/disk/bob/gpm/n1_oly_ocean_ku/classify/ex_data/gpm-ku_cfad_total.txt'
outFileKuRebin = '/home/disk/bob/gpm/n1_oly_ocean_ku/classify/ex_data/gpm-ku-rebin_cfad_total.txt'
rawDataKa = np.loadtxt(inFileKa)
rawDataKu = np.loadtxt(inFileKu)

#redimension rawDataKu - cut levels in half (from 176 to 88) to match dims of rawDataKa
rows = rawDataKa.shape[0]
rawDataKuRebin = np.empty(rawDataKa.shape,dtype=float)
for i in range(0,rows):
    rawDataKuRebin[i,:] = rawDataKu[2*i,:] + rawDataKu[2*i+1,:]
#fid = open(outFileKuRebin,'w')
#for irow in range(0,rows):
#    rawDataKuRebin[irow,:].tofile(fid,sep="\t",format="%d")
#    fid.write("\n")
#fid.close()

#normalize data
rawDataKaNorm = rawDataKa / np.max(rawDataKa)
rawDataKuNorm = rawDataKuRebin / np.max(rawDataKuRebin)

#diff the two input arrays
rawDataNorm = rawDataKaNorm - rawDataKuNorm

#bin dimensions
xcenter = np.linspace(0.5,49.5,50)
ycenter = np.linspace(0.0,21.875,88)
X_grid, Y_grid = np.meshgrid(xcenter,ycenter)

#define colorbar range
#cbar_levels = np.linspace(0.05,1,20)
#cbar_levels = np.linspace(-0.6,0.3,14)
cbar_levels = np.linspace(-0.4,0.4,17)

#plot cfad
fig = plt.figure()
fig.set_size_inches(6,6)
ax0 = fig.add_subplot(1,1,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title('GPM Ka-Ku Normalized CFAD (DBZ) - Over Ocean')
cfad = ax0.contourf(X_grid, Y_grid, rawDataNorm, levels=cbar_levels,cmap='jet')
#cfad = ax0.contourf(X_grid, Y_grid, rawDataNorm, levels=cbar_levels,cmap='Spectral')

ax0.set_xlabel('Reflectivity (dBZ)')
ax0.set_ylabel('Height (km)')
ax0.set_ylim(0,8)
ax0.set_xlim(0,50)
ax0.grid()

#cbar = plt.colorbar(cfad,cax,ticks=np.linspace(0.05,1,20))
cbar = plt.colorbar(cfad,cax,ticks=np.linspace(-0.4,0.4,17))
cbar.ax.set_ylabel(r'Frequency')

plt.savefig(outdir+'gpm_cfad_ka-ku.png', bbox_inches='tight')

