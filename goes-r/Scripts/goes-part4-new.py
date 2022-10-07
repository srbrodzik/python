#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:49:18 2019

@author: brodzik
"""

# GNC-A Blog GOES-16 Python Tutorial: Part IV - Adding a background, changing transparency and adding a land/sea maak
 
# Required libraries
import matplotlib.pyplot as plt # Import the Matplotlib package
#%matplotlib inline
from netCDF4 import Dataset # Import the NetCDF Python interface
from mpl_toolkits.basemap import Basemap # Import the Basemap toolkit
import numpy as np # Import the Numpy package
import os
from convertRadiance2 import radianceToRefl
from convertRadiance2 import radianceToTb

# Path to the GOES-R image file
#path = '/home/disk/shear2/brodzik/python/goes-r/DATA/L2/FD/Ch13/OR_ABI-L2-CMIPF-M3C13_G16_s20190600000311_e20190600011089_c20190600011164.nc'
path = '/home/disk/shear2/brodzik/python/goes-r/DATA/L1b/FD/Ch13/OR_ABI-L1b-RadF-M3C13_G16_s20190600000311_e20190600011089_c20190600011142.nc'
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
        kappa = nc.variables['kappa0'][:]
        data = radianceToRefl(rad,kappa)
    else:
        fk1 = nc.variables['planck_fk1'][:]
        fk2 = nc.variables['planck_fk2'][:]
        bc1 = nc.variables['planck_bc1'][:]
        bc2 = nc.variables['planck_bc2'][:]
        data = radianceToTb(rad,fk1,fk2,bc1,bc2)
elif Level == '2':
    data = nc.variables['CMI'][:]

# Create the basemap reference for the Satellite Projection
# resolution options: ‘c’ (crude), ‘l’ (low), ‘i’ (intermediate), ‘h’ (high) and ‘f’ (full).
#bmap = Basemap(projection='geos', lon_0=-89.5, lat_0=0.0, satellite_height=35786023.0, ellps='GRS80', resolution = 'i')
bmap = Basemap(projection='geos', lon_0=-75.2, lat_0=0.0, satellite_height=35786023.0, ellps='GRS80')
#
# Add a background
#bmap.bluemarble()
#bmap.shadedrelief()
bmap.etopo()

# Add an ocean mask - useful for land-only products
#bmap.drawlsmask(ocean_color='aqua',lakes=True)

# Add a land mask - useful for ocean-only products
bmap.drawlsmask(land_color='chartreuse')

# Plot GOES-16 channel with transparency
bmap.imshow(data, origin='upper', vmin=170, vmax=378, cmap='Greys', alpha=0.7)

# Plot GOES-16 Channel using 170 and 378 as the temperature thresholds with no background
#bmap.imshow(data, origin='upper', vmin=170, vmax=378, cmap='Greys')
#bmap.imshow(data, origin='upper', vmin=170, vmax=320, cmap='jet')
#bmap.imshow(data, origin='upper', vmin=170, vmax=360, cmap='hsv')

# Draw the coastlines, countries, parallels and meridians
bmap.drawcoastlines(linewidth=0.5, linestyle='solid', color='black')
bmap.drawcountries(linewidth=0.5, linestyle='solid', color='black')
bmap.drawparallels(np.arange(-90.0, 90.0, 10.0), linewidth=0.5, color='white')
bmap.drawmeridians(np.arange(0.0, 360.0, 10.0), linewidth=0.5, color='white')
 
# Insert the legend at the bottom
#bmap.colorbar(location='bottom', label='Brightness Temperature [K]')
 
# Add a title to the plot
#plt.title("GOES-16 ABI Band " + Band + "\n Scan from " + Start_Formatted + "\n to " + End_Formatted)
#plt.title("GOES-16 ABI Band " + Band + " - Blue Marble Background\n Scan from " + Start_Formatted + "\n to " + End_Formatted)
plt.title("GOES-16 ABI Band " + Band + " - Shaded Relief Background\n Scan from " + Start_Formatted + "\n to " + End_Formatted, fontsize=8)
#plt.title("GOES-16 ABI Band " + Band + " - Etopo Background\n Scan from " + Start_Formatted + "\n to " + End_Formatted)

# Read some variables from the NetCDF header in order to use it in the plot
geo_extent = nc.variables['geospatial_lat_lon_extent']
 
center = str(geo_extent.geospatial_lon_center)
west = str(geo_extent.geospatial_westbound_longitude)
east = str(geo_extent.geospatial_eastbound_longitude)
north = str(geo_extent.geospatial_northbound_latitude)
south = str(geo_extent.geospatial_southbound_latitude)
 
# Put the information retrieved from the header in the final image
#plt.text(-300000,300000,'Geospatial Extent \n' + west + '°W \n' + east + '°E \n' + north + '°N \n' + south + '°S \n' + 'Center = ' + center + '°', fontsize = 7)
plt.text(-300000,300000,'Geospatial Extent \n' + west + 'W \n' + east + 'E \n' + north + 'N \n' + south + 'S \n' + 'Center = ' + center , fontsize = 5)

# Export result
DPI = 150
#plt.savefig('/home/disk/shear2/brodzik/python/goes-r/PLOTS/goes-part4_Ch13_landMask0.7.png', dpi=DPI, bbox_inches='tight', pad_inches=0)

# Show the plot
plt.show()
