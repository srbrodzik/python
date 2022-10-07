#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#======================================================================================================
# GNC-A Blog Python Tutorial: Part VI
#======================================================================================================# Required libraries ==================================================================================
import sys

import matplotlib.pyplot as plt # Import the Matplotlib package
from mpl_toolkits.basemap import Basemap # Import the Basemap toolkit&amp;amp;amp;amp;amp;amp;amp;amp;amp;lt;/pre&amp;amp;amp;amp;amp;amp;amp;amp;amp;gt;
import numpy as np # Import the Numpy package
 
from remap import remap # Import the Remap function
 
from cpt_convert import loadCPT # Import the CPT convert function
from matplotlib.colors import LinearSegmentedColormap # Linear interpolation for color maps

import datetime # Library to convert julian day to dd-mm-yyyy

from miscFunctions import centerWavelength
#======================================================================================================

# Read input arguments - if you use this, comment out line 31 below
#if len(sys.argv) != 2:
#    print >>sys.stderr, "Useage: ",sys.argv[0]," [full path for input file]"
#    quit()
#path = sys.argv[1]

 
# Load the Data =======================================================================================
# Path to the GOES-16 image file
path = '/home/disk/shear2/brodzik/python/goes-r/DATA/OR_ABI-L2-CMIPF-M3C13_G16_s20190600000311_e20190600011089_c20190600011164.nc'
 
# Choose the visualization extent (min lon, min lat, max lon, max lat)
# Amazon Region
#extent = [-90.0, -40.0, -20.0, 10.0]
# North & South America
#extent = [-160.0, -80.0, 10.0, 80.0]
# Columbia
#extent = [-85.0, -5.0, -60.0, 12.0]
# USA
extent = [-130.0, 20.0, -60.0, 55.0]

# Choose the image resolution (the higher the number the faster the processing is)
resolution = 2.0
 
# Call the reprojection funcion
grid = remap(path, extent, resolution, 'HDF5')
 
# Read the data returned by the function
data = grid.ReadAsArray()
#======================================================================================================
 
# Plot the Data =======================================================================================
# Create the basemap reference for the Rectangular Projection
bmap = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)
 
# Draw the countries and Brazilian states shapefiles
bmap.readshapefile('/home/disk/shear2/brodzik/python/goes-r/DATA/ShapeFiles/Brazil/BRA_adm1','BRA_adm1',linewidth=0.50,color='darkslategray')
bmap.readshapefile('/home/disk/shear2/brodzik/python/goes-r/DATA/ShapeFiles/Countries_Shape/ne_10m_admin_0_countries','ne_10m_admin_0_countries',linewidth=0.50,color='darkslategray')
 
# Draw parallels and meridians
bmap.drawparallels(np.arange(-90.0, 90.0, 5.0), linewidth=0.25, color='white', labels=[True,False,False,True])
bmap.drawmeridians(np.arange(0.0, 360.0, 5.0), linewidth=0.25, color='white', labels=[True,False,False,True])
 
# Converts a CPT file to be used in Python
cpt = loadCPT('/home/disk/shear2/brodzik/python/colors/IR4AVHRR6.cpt')
# Makes a linear interpolation
cpt_convert = LinearSegmentedColormap('cpt', cpt)
 
# Plot the GOES-16 channel with the converted CPT colors
bmap.imshow(data, origin='upper', cmap=cpt_convert, vmin=170, vmax=378)
 
# Search for the GOES-16 channel in the file name
Band = (path[path.find("M3C")+3:path.find("_G16")])
# Search for the Scan start in the file name
Start = (path[path.find("_s")+2:path.find("_e")])
Start_Formatted = Start[0:4] + " Day " + Start[4:7] + " - " + Start [7:9] + ":" + Start [9:11] + ":" + Start [11:13] + "." + Start [13:14] + " UTC"
# Search for the Scan end in the file name
End = (path[path.find("_e")+2:path.find("_c")])
End_Formatted = End[0:4] + " Day " + End[4:7] + " - " + End [7:9] + ":" + End [9:11] + ":" + End [11:13] + "." + End [13:14] + " UTC"

# Convert julian days to DD-MMM-YYYY
year = int(Start[0:4])
dayjulian = int(Start[4:7]) - 1
daynormal = datetime.datetime(year,1,1) + datetime.timedelta(dayjulian)
daynormal2 = daynormal.strftime('%A %d %b %Y')
time = Start [7:9] + ":" + Start [9:11] + ":" + Start [11:13] + "." + Start [13:14] + " UTC"
 
# Add a title to the plot
#plt.title("GOES-16 ABI Band " + Band + "\n Scan from " + Start_Formatted + " to " + End_Formatted, fontsize=8)
#plt.title("GOES-16 ABI Band " + Band + " " + Center_WL + " " + daynormal2 + " " + time)
plt.title("GOES-16 ABI Band " + Band + " " + centerWavelength(Band) + " " + daynormal2 + " " + time)

# Insert the legend at the bottom
bmap.colorbar(location='right', label='Brightness Temperature [K]')
 
# Show the plot
plt.show()
#======================================================================================================;