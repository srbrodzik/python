'''
plot_mrr_cfad.py
Make CFAD of mrr data from Atmos Rivers
Not sure what's going to come up...
for BB:
-Add 86.7 m to everything for CRN elevation
-Otherwise keep original
'''
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pyart
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cPickle as pickle
import pdb
import numpy as np


infile = '/home/disk/meso-home/jzagrod/Olympex/MRR/CFAD/Test/Data/fishery_atmos_rivers_warm_sector.p'
outdir = '/home/disk/meso-home/jzagrod/Olympex/MRR/CFAD/Test/Plot/'
cfad = pickle.load( open(infile, "rb" ) )

time = cfad[0]
ranges = cfad[1]
dbz = cfad[2]
vel = cfad[3]
bb_height = cfad[4]
bb_int = int(np.mean(bb_height[np.where((bb_height > 0))[0]])+86.7)

#Correct ranges for BB:
mean_bb_h = np.mean(bb_height[np.where((bb_height > 0))[0]])
for b,bb_h in enumerate(bb_height):
    ranges[b,:] += 86.7
    if bb_h > 0:
        ranges[b,:] += mean_bb_h - bb_h


xedges_dbz = np.linspace(10,50,41)
xcenter_dbz = np.linspace(10.5,49.5,40)
xedges_vel = np.linspace(0,10,51)
xcenter_vel = np.linspace(0.1,9.9,50)
yedges = np.linspace(50,3150,32)
ycenter = np.linspace(100,3100,31)

#Find mean profiles
mean_dbz = np.zeros(31)
mean_vel = np.zeros(31)
dbzall = np.ravel(dbz)
rangesall = np.ravel(ranges)
for i,level in enumerate(ycenter):
    mean_dbz[i] = np.mean(dbzall[np.where((rangesall >= yedges[i]) & (rangesall < yedges[i+1]))[0]])

dbz_cfad, xedges_dbz, yedges = np.histogram2d(np.ravel(dbz), np.ravel(ranges), bins=(xedges_dbz, yedges))
vel_cfad, xedges_vel, yedges = np.histogram2d(np.ravel(vel), np.ravel(ranges), bins=(xedges_vel, yedges))

dbz_final = np.zeros(np.shape(dbz_cfad))
vel_final = np.zeros(np.shape(vel_cfad))

for i,cf in enumerate(dbz_cfad[0,:]):
    hist = dbz_cfad[:,i]  
    dbz_final[:,i] = hist/np.sum(hist)*100.
    hist2 = vel_cfad[:,i]  
    vel_final[:,i] = hist2/np.sum(hist2)*100.
    

#plot CFAD...2-panel
X_dbz, Y_dbz = np.meshgrid(xcenter_dbz, ycenter)
dbz_final = np.rot90(dbz_final)
dbz_final = np.flipud(dbz_final)

X_vel, Y_vel = np.meshgrid(xcenter_vel, ycenter)
vel_final = np.rot90(vel_final)
vel_final = np.flipud(vel_final)

dbzlevels = np.linspace(0.25,15,15*4)
vellevels = np.linspace(0.25,35,35*4)

labelname = '# obs (total = '+str(len(time))+')'
fig = plt.figure()
fig.set_size_inches(12,6)
ax0 = fig.add_subplot(1,2,1) 
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax0.set_title('Fishery MRR ARs (Warm Sector)')
histplot = ax0.contourf(X_dbz, Y_dbz, dbz_final, label = labelname,levels=dbzlevels,cmap='jet')
#ax0.plot(mean_dbz,ycenter, 'k-',label="Mean",color="black",lw=3)
ax0.plot(xedges_dbz,np.zeros(np.shape(xedges_dbz))+bb_int, 'k--',label="Mean",color="black",lw=2)

ax0.set_xlabel('Reflectivity (dBZ)')
ax0.set_ylabel('Normalized height (m)')
ax0.set_ylim(0,3100)
ax0.set_xlim(10,45)
ax0.grid()

ax0.text(0.97, 0.95,'# 1-min obs: '+str(len(time)), horizontalalignment='right',
      verticalalignment='center',
      transform=ax0.transAxes)
cbar = plt.colorbar(histplot,cax,ticks=np.linspace(1,15,15))
cbar.ax.set_ylabel(r'Frequency (%)')

ax1 = fig.add_subplot(1,2,2) 
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="3%", pad=0.15)
ax1.set_title(' ')
histplot = ax1.contourf(X_vel, Y_vel, vel_final, label = labelname,levels=vellevels,cmap='jet')
#ax1.plot(xedges_vel,np.zeros(np.shape(xedges_vel))+bb_int, 'k--',label="Mean",color="black",lw=2)
ax1.plot(mean_vel,ranges[0,:], 'k-',label="Mean",color="black",lw=3)

ax1.set_xlabel('Velocity (m s$^{-1}$)')
ax1.set_ylabel('Normalized height (m)')
ax1.set_ylim(0,3100)
ax1.set_xlim(0,10)
ax1.grid()
'''
ax1.text(0.07, 0.95,'# 1-min obs: '+str(len(time)), horizontalalignment='right',
      verticalalignment='center',
      transform=ax0.transAxes)
'''
cbar = plt.colorbar(histplot,cax,ticks=np.linspace(2,36,18))
cbar.ax.set_ylabel(r'Frequency (%)')



plt.tight_layout()

plt.savefig(outdir+'cfad_fishery_atmos_rivers_warm_sector.png', bbox_inches='tight')



pdb.set_trace()
