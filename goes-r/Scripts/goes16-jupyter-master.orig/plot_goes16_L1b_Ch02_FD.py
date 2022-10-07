#%matplotlib
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from pyproj import Proj

#Read in level 1b GOES-16 NetCDF
g16nc = Dataset('/home/disk/data/satellite/GOES/GRB16/ABI/CONUS/Channel02/20190430/OR_ABI-L1b-RadC-M6C02_G16_s20191201746225_e20191201748598_c20191201749039.nc', 'r')
radiance = g16nc.variables['Rad'][:]

#Convert Radiance to Reflectance
#It is useful to convert the radiance values from an absolute scale into a relative scale to make them
#easier to deal with in the following steps. We will use the example formula published by NOAA at
#http://www.goes-r.gov/products/ATBDs/baseline/Imagery_v2.0_no_color.pdf on page 21.

# Define some constants needed for the conversion. From the pdf linked above
Esun_Ch_02 = 663.274497
d2 = 0.3

# Apply the formula to convert radiance to reflectance
ref = (radiance * np.pi * d2) / Esun_Ch_02

# Make sure all data is in the valid data range
ref = np.maximum(ref, 0.0)
ref = np.minimum(ref, 1.0)

#Gamma correction
#You'll notice the image above looks very dark, that is because the values are in linear units.
#We do a simple gamma correction to adjust this and brighten the image. This makes it easier to see
#features in the data.

# Apply the formula to adjust reflectance gamma
ref_gamma = np.sqrt(ref)

# Plot gamma adjusted reflectance
#fig = plt.figure(figsize=(6,6),dpi=200)
#im = plt.imshow(ref_gamma, vmin=0.0, vmax=1.0, cmap='Greys_r')
#cb = fig.colorbar(im, orientation='horizontal')
#cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
#cb.set_label('Reflectance Gamma Adjusted')
#plt.show()

#Georeferenced GOES-16 Image

# Height of the satellite's orbit
height = g16nc.variables['goes_imager_projection'].perspective_point_height

# Longitude of the satellite's orbit
lon = g16nc.variables['goes_imager_projection'].longitude_of_projection_origin

# Which direction do we sweep in?
sweep = g16nc.variables['goes_imager_projection'].sweep_angle_axis

X = g16nc.variables['x'][:] * height
Y = g16nc.variables['y'][:] * height

#Lets use Basemap to map state, coastline & country borders now.
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(ref_gamma, vmin=0.0, vmax=1.0, cmap='Greys_r')
m = Basemap(resolution='l', projection='geos', lon_0=lon,
            llcrnrx=X.min(), llcrnry=Y.min(), urcrnrx=X.max(), urcrnry=Y.max())
m.imshow(np.flipud(ref_gamma), vmin=0.0, vmax=1.0, cmap='Greys_r')
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.drawparallels(np.arange(10.0, 60.0, 10.0), linewidth=0.5, color='blue')
m.drawmeridians(np.arange(230.0, 310.0, 10.0), linewidth=0.5, color='blue')
cb = fig.colorbar(im, orientation='horizontal')
cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('Reflectance Gamma Adjusted')
plt.show()

g16nc.close()
g16nc = None

