#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:14:22 2019

@author: brodzik
"""

#GNC-A Blog GOES-16 Python Tutorial: Part III - Getting Information from the file name and header

# Required libraries
import matplotlib.pyplot as plt # Import the Matplotlib package
from netCDF4 import Dataset # Import the NetCDF Python interface
from mpl_toolkits.basemap import Basemap # Import the Basemap toolkit
import numpy as np # Import the Numpy package
import os
from convertRadiance import radianceToRefl
from convertRadiance import radianceToTb

# Path to the GOES-R image file
path = '/home/disk/shear2/brodzik/python/goes-r/DATA/OR_ABI-L2-CMIPF-M3C13_G16_s20190600000311_e20190600011089_c20190600011164.nc'
#path = '/home/disk/shear2/brodzik/python/goes-r/DATA/OR_ABI-L1b-RadF-M3C13_G16_s20190600000311_e20190600011089_c20190600011142.nc'
fname = os.path.basename(path)

# Filenaming convention: 
#    OR_ABI-L2_CMIP[C/F/M]-M[#][ch]_G[satno]_s[YYYYJJJHHMMSSs]_e[YYYYJJJHHMMSSs]_c[YYYYJJJHHMMSSs].nc
#    [C/F/M] - CONUS/FULLDISK/MESOSCALE
#    [#] - ABI Mode Number
#          Mode 3  15 min Full Disk, 5 min CONUS, 1 min/30 sec Mesoscale(s)
#          Mode 4   5 min Full Disk
#          Mode 6  10 min Full Disk, 5 min CONUS, 1 min/30 sec Mesoscale(s)
#    [ch] - Channel number from '01' to '16'
#    [satno] - 16 (or 17)
#    obs_start/obs_end/file_creation dates

# Search for the file level in the file name
Level = (fname[fname.find("ABI-L")+5])
# Search for the GOES-R channel in the file name
#Band = (fname[fname.find("M3C")+3:fname.find("_G16")])
Band = (fname[fname.find("_G16")-2:fname.find("_G16")])
# Search for the Scan Start in the file name
Start = (fname[fname.find("s")+1:fname.find("_e")])
Start_Formatted = Start[0:4] + " Day " + Start[4:7] + " - " + Start [7:9] + ":" + Start [9:11] + ":" + Start [11:13] + "." + Start [13:14] + " UTC"
# Search for the Scan End in the file name
End = (fname[fname.find("e")+1:fname.find("_c")])
End_Formatted = End[0:4] + " Day " + End[4:7] + " - " + End [7:9] + ":" + End [9:11] + ":" + End [11:13] + "." + End [13:14] + " UTC"

# Open the file using the NetCDF4 library
nc = Dataset(path)

# Extract the Brightness Temperature or radiance values from the NetCDF
if Level == '1':
    rad = nc.variables['Rad'][:] 
    if Band <= '06':
        data = radianceToRefl(rad,Band)
    else:
        data = radianceToTb(rad,Band)
elif Level == '2':
    data = nc.variables['CMI'][:]

#Suggested update from tutorial
#If these lines are used, comment out the DPI command at bottom of code
DPI = 250
ax = plt.figure(figsize=(2000/float(DPI), 2000/float(DPI)), frameon=True, dpi=DPI)

# Create the basemap reference for the Satellite Projection
# resolution options: ‘c’ (crude), ‘l’ (low), ‘i’ (intermediate), ‘h’ (high) and ‘f’ (full).
#bmap = Basemap(projection='geos', lon_0=-75.2, lat_0=0.0, satellite_height=35786023.0, ellps='GRS80', resolution = 'i')
bmap = Basemap(projection='geos', lon_0=-75.2, lat_0=0.0, satellite_height=35786023.0, ellps='GRS80')
#
# Plot GOES-16 Channel using 170 and 378 as the temperature thresholds
#bmap.imshow(data, origin='upper', vmin=170, vmax=378, cmap='Greys')
#map.imshow(data, origin='upper', vmin=170, vmax=320, cmap='jet')
bmap.imshow(data, origin='upper', vmin=170, vmax=360, cmap='hsv')

# Draw the coastlines, countries, parallels and meridians
bmap.drawcoastlines(linewidth=0.5, linestyle='solid', color='black')
bmap.drawcountries(linewidth=0.5, linestyle='solid', color='black')
bmap.drawparallels(np.arange(-90.0, 90.0, 10.0), linewidth=0.5, color='white')
bmap.drawmeridians(np.arange(0.0, 360.0, 10.0), linewidth=0.5, color='white')
 
# Insert the legend at the bottom
bmap.colorbar(location='bottom', label='Brightness Temperature [K]')
 
# Add a title to the plot
plt.title("GOES-16 ABI Band " + Band + "\n Scan from " + Start_Formatted + "\n to " + End_Formatted)

# Read some variables from the NetCDF header in order to use it in the plot
geo_extent = nc.variables['geospatial_lat_lon_extent']
 
center = str(geo_extent.geospatial_lon_center)
west = str(geo_extent.geospatial_westbound_longitude)
east = str(geo_extent.geospatial_eastbound_longitude)
north = str(geo_extent.geospatial_northbound_latitude)
south = str(geo_extent.geospatial_southbound_latitude)
 
# Put the information retrieved from the header in the final image
#plt.text(-300000,300000,'Geospatial Extent \n' + west + '°W \n' + east + '°E \n' + north + '°N \n' + south + '°S \n' + 'Center = ' + center + '°', fontsize = 7)
plt.text(-300000,300000,'Geospatial Extent \n' + west + 'W \n' + east + 'E \n' + north + 'N \n' + south + 'S \n' + 'Center = ' + center , fontsize = 7)

# Export result
#DPI = 150
plt.savefig('/home/disk/shear2/brodzik/python/goes-r/PLOTS/goes-part3-hsv-fromL2.png', dpi=DPI, bbox_inches='tight', pad_inches=0)

# Show the plot
plt.show()
