#!/usr/bin/python3

# Based on example here:
# https://unidata.github.io/MetPy/latest/examples/cross_section.html#sphx-glr-examples-cross-section-py

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.ticker import *
import matplotlib.ticker as ticker
import numpy as np
import xarray as xr

import metpy.calc as mpcalc
from metpy.interpolate import cross_section
from metpy.plots import ctables

infile = '/home/disk/bob/impacts/raw/gpm_ku/202002/GPM2Ku6_uw3_20200201.143222_to_20200201.143921_033679_IMP.nc'
verbose = False

# Read gridded data and squeeze down size-one time dimension
data = xr.open_dataset(infile)
data = data.metpy.parse_cf().squeeze()
if verbose:
    print(data)

# Define start and end points
start = (35.341034, -73.51219)
end = (37.64104, -71.7499)

# Get the cross-section
cross = cross_section(data, start, end, steps=1000, interp_type='linear').set_coords(('lat', 'lon'))
if verbose:
    print(cross)

# Define the figure object and primary axes
fig = plt.figure(1, figsize=(10., 6.))
ax = plt.axes()

# Plot refl using contourf
refl_contour = ax.contourf(cross['lat'], cross['alt'], cross['refl'],
                           levels=np.arange(-10, 50.01, 0.5),
                           cmap=ctables.registry.get_colortable('NWSStormClearReflectivity'))
refl_colorbar = fig.colorbar(refl_contour, ticks=MaxNLocator(7))

# Tailor y-axis to match other plot
ax.set_ylim(0,12)
tick_spacing = 0.5
ax.yaxis.set_minor_locator(ticker.MultipleLocator(tick_spacing))

# Set axes and colorbar labels
ax.set_ylabel('Altitude (km ASL)')
ax.set_xlabel('Latitude (degrees north)')
refl_colorbar.set_label('GPM Ku-band (dBZ)')
ax.grid()

plt.show()
