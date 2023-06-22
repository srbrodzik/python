# Brody Fuchs, CSU, Nov 2017

# this will take in rays and plot them from the ascope program

import numpy as np 
import glob
import os
import matplotlib.pyplot as plt 

base_path = '/home/rmet/SPURS2/ascope_data'

fname_ppi = '%s/ppi_ray.npz' % (base_path)
fname_rhi = '%s/rhi_ray.npz' % (base_path)

ppi_ray = np.load(fname_ppi)
rhi_ray = np.load(fname_rhi)

rhi_range_corr = rhi_ray['y'] + 20.0*np.log10(rhi_ray['x'])

rrc_diff = rhi_range_corr-ppi_ray['y'][:len(rhi_range_corr)]

print np.nanmean(rrc_diff), np.nanstd(rrc_diff)

fig, ax = plt.subplots(1,1, figsize=(8,6))

ax.plot(ppi_ray['x'], ppi_ray['y'], label='ppi', linewidth=0.7, color='purple')
ax.plot(rhi_ray['x'], rhi_ray['y'], label='rhi', linewidth=0.7, color='black')

ax.plot(rhi_ray['x'], rhi_range_corr, label='rhi corrected', linewidth=0.7, color='gray')
#ax.plot()

ax2 = ax.twinx()

ax2.plot(rhi_ray['x'], rrc_diff, color='red', alpha=0.5, linewidth=0.5)

ax2.set_ylim(0, 50)

ax2.set_ylabel('dBZ diff', color='red')

ax.grid(True)
ax.set_xlabel('Range (km)')
ax.set_ylabel('dBZ')

ax.legend(loc='best')

plt.tight_layout()

plt.savefig('ray_comparison_test.png')