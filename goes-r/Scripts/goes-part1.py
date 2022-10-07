# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# GNC-A Blog GOES-16 Python Tutorial: Part I - Extracting Values and Basic Plotting
 
# Required libraries
import matplotlib.pyplot as plt
from netCDF4 import Dataset
 
# Path to the GOES-R image file
path = '/home/disk/shear2/brodzik/python/goes-r/DATA/OR_ABI-L2-CMIPF-M3C13_G16_s20190600000311_e20190600011089_c20190600011164.nc'
# Try using L1b data
#path = '/home/disk/shear2/brodzik/python/goes-r/DATA/OR_ABI-L1b-RadF-M3C13_G16_s20190600000311_e20190600011089_c20190600011142.nc'

# Open the file using the NetCDF4 library
nc = Dataset(path)
 
# L2: Extract the Brightness Temperature values from the NetCDF
data = nc.variables['CMI'][:]
# L1b: Extract the Radiance values from the NetCDF
#data = nc.variables['Rad'][:]

# Show data
plt.imshow(data, cmap='Greys')
#plt.imshow(T, cmap='Greys')