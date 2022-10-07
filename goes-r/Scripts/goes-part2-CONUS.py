#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 13:45:47 2019

@author: brodzik
"""

# GNC-A Blog GOES-16 Python Tutorial: Part II
 
# Required libraries
import matplotlib.pyplot as plt # Import the Matplotlib package
from netCDF4 import Dataset # Import the NetCDF Python interface
from mpl_toolkits.basemap import Basemap # Import the Basemap toolkit
import numpy as np # Import the Numpy package
 
# Path to the GOES-R image file
path = '/home/disk/shear2/brodzik/python/goes-r/DATA/OR_ABI-L2-CMIPC-M6C13_G16_s20191000016194_e20191000018579_c20191000019022.nc'
 
# Open the file using the NetCDF4 library
nc = Dataset(path)
 
# Extract the Brightness Temperature values from the NetCDF
data = nc.variables['CMI'][:] 
 
#Suggested update from tutorial
#If these lines are used, comment out the DPI command at bottom of code
DPI = 150
ax = plt.figure(figsize=(2000/float(DPI), 2000/float(DPI)), frameon=True, dpi=DPI)

# Create the basemap reference for the Satellite Projection
# resolution options: ‘c’ (crude), ‘l’ (low), ‘i’ (intermediate), ‘h’ (high) and ‘f’ (full).
#bmap = Basemap(projection='geos', lon_0=-75.2, lat_0=0.0, satellite_height=35786023.0, ellps='GRS80', resolution = 'i')
#bmap = Basemap(projection='geos', lon_0=-75.2, lat_0=0.0, satellite_height=35786023.0, ellps='GRS80')
bmap = Basemap(projection='merc', llcrnrlat=10.0, urcrnrlat=60.0, llcrnrlon=-155, urcrnrlon=-50, lat_ts=5, resolution='c')

# Plot GOES-16 Channel using 170 and 378 as the temperature thresholds
#bmap.imshow(data, origin='upper', vmin=170, vmax=378, cmap='Greys')
#bmap.imshow(data, origin='upper', vmin=170, vmax=325, cmap='jet')
bmap.imshow(data, origin='upper', vmin=170, vmax=360, cmap='hsv')
 
# Draw the coastlines, countries, parallels and meridians
#bmap.drawcoastlines(linewidth=0.3, linestyle='solid', color='black')
#bmap.drawcountries(linewidth=0.3, linestyle='solid', color='black')
#bmap.drawparallels(np.arange(-90.0, 90.0, 10.0), linewidth=0.1, color='white')
#bmap.drawmeridians(np.arange(0.0, 360.0, 10.0), linewidth=0.1, color='white')
bmap.drawcoastlines(linewidth=0.5, linestyle='solid', color='black')
bmap.drawcountries(linewidth=0.5, linestyle='solid', color='black')
bmap.drawparallels(np.arange(10.0, 60.0, 5.0), linewidth=0.3, color='white')
bmap.drawmeridians(np.arange(-155.0, -50.0, 5.0), linewidth=0.3, color='white')

# Insert the legend
bmap.colorbar(location='bottom', label='Brightness Temperature [K]')
 
# Export result
#DPI = 300   # orig
#DPI = 150   # smaller
plt.savefig('/home/disk/shear2/brodzik/python/goes-r/PLOTS/goes-part2-hsv-CONUS-merc.png', dpi=DPI, bbox_inches='tight', pad_inches=0)

# Show the plot
plt.show()