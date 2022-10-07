%matplotlib
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#Read in level 1b GOES-16 NetCDF
g16nc = Dataset('OR_ABI-L1b-RadM1-M3C02_G16_s20171931811268_e20171931811326_c20171931811356.nc', 'r')
radiance = g16nc.variables['Rad'][:]
g16nc.close()
g16nc = None

#Make initial Radiance plot
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(radiance, cmap='Greys_r')
cb = fig.colorbar(im, orientation='horizontal')
cb.set_ticks([1, 100, 200, 300, 400, 500, 600])
cb.set_label('Radiance (W m-2 sr-1 um-1)')
plt.show()

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

# Plot reflectance
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(ref, vmin=0.0, vmax=1.0, cmap='Greys_r')
cb = fig.colorbar(im, orientation='horizontal')
cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('Reflectance')
plt.show()

#Reading in Level 2 GOES NetCDF data
#So the Level 1b data contains radiance values, and we just demonstrated the process to convert this
#data to reflectance. However, lets say you just want to get straight to the good stuff and not go
#through this process? Thankfully, NOAA provides the reflectance & brightness temperatures precomputed
#as the Level 2 product. We will work with those from here on out.

# Read in Level 2 NetCDF
g16nc = Dataset('OR_ABI-L2-CMIPM1-M3C02_G16_s20171931811268_e20171931811326_c20171931811393.nc', 'r')
ref_ch2 = g16nc.variables['CMI'][:]
g16nc.close()
g16nc = None

# Plot reflectance
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(ref_ch2, vmin=0.0, vmax=1.0, cmap='Greys_r')
cb = fig.colorbar(im, orientation='horizontal')
cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('Level 2 Reflectance')
plt.show()

#Gamma correction
#You'll notice the image above looks very dark, that is because the values are in linear units.
#We do a simple gamma correction to adjust this and brighten the image. This makes it easier to see
#features in the data.

# Apply the formula to adjust reflectance gamma
ref_gamma = np.sqrt(ref_ch2)

# Plot gamma adjusted reflectance
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(ref_gamma, vmin=0.0, vmax=1.0, cmap='Greys_r')
cb = fig.colorbar(im, orientation='horizontal')
cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('Reflectance Gamma Adjusted')
plt.show()

#Pseudo-True Color Image
#GOES-16 Contains a Blue band and a Red band but no Green band. Fortunately, we can derive a linear
#relationship between the blue, red, and veggie (near IR) bands to approximate a green band. We will
#do that here.

# Load Channel 1 - Blue Visible & do gamma adjustment
g16nc = Dataset('OR_ABI-L2-CMIPM1-M3C01_G16_s20171931811268_e20171931811326_c20171931811382.nc', 'r')
ref_1 = g16nc.variables['CMI'][:]
g16nc.close()
g16nc = None
ref_gamma_1 = np.sqrt(ref_1)

# Plot gamma adjusted reflectance channel 1
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(ref_gamma_1, vmin=0.0, vmax=1.0, cmap='Greys_r')
cb = fig.colorbar(im, orientation='horizontal')
cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('Ch01 - Reflectance')
plt.show()

# Load Channel 3 - Veggie Near IR & do gamma adjustment
g16nc = Dataset('OR_ABI-L2-CMIPM1-M3C03_G16_s20171931811268_e20171931811326_c20171931811389.nc', 'r')
ref_3 = g16nc.variables['CMI'][:]
g16nc.close()
g16nc = None
ref_gamma_3 = np.sqrt(ref_3)

# Plot gamma adjusted reflectance channel 3
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(ref_gamma_3, vmin=0.0, vmax=1.0, cmap='Greys_r')
cb = fig.colorbar(im, orientation='horizontal')
cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('Ch03 - Reflectance')
plt.show()

#On GOES-16, Band 2 (Red Visible) is 500-meter resolution while Band 1 (Blue Visible) and Band 3
#(Veggie IR) are 1,000-meter resolution. In order to combine the 3 bands into an RGB image we
#will first resample bands 2 to 1000-meter resolution.

#Rebin function from
#https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

ref_gamma_2 = rebin(ref_gamma, [1000, 1000])

#Lets make a geocolor image using the veggie near IR band directly in place of green and see how that looks.
geocolor = np.stack([ref_gamma_2, ref_gamma_3, ref_gamma_1], axis=2)
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(geocolor)
plt.title('GeoColor - Red - Veggie - Blue')
plt.show()

#As you can see, the green channel overpowers everything else. We will use a simple linear
#relationship to adjust the green values to create a psuedo-green channel and produce a better
#true color image.

#Derived from Planet Labs data, CC > 0.9
ref_gamma_true_green = 0.48358168 * ref_gamma_2 + 0.45706946 * ref_gamma_1 + 0.06038137 * ref_gamma_3

truecolor = np.stack([ref_gamma_2, ref_gamma_true_green, ref_gamma_1], axis=2)
fig = plt.figure(figsize=(6,6),dpi=200)
im = plt.imshow(truecolor)
plt.title('TrueColor - Red - Psuedo-Green - Blue')
plt.show()

#Much better!  Now lets georeference the data so we can plot it with Basemap

#Georeferenced GOES-16 Image
from mpl_toolkits.basemap import Basemap
from pyproj import Proj

#We need to load in some more variables from the netCDF file so we will reopen it to grab them now.
g16nc = Dataset('OR_ABI-L2-CMIPM1-M3C03_G16_s20171931811268_e20171931811326_c20171931811389.nc', 'r')

# Height of the satellite's orbit
height = g16nc.variables['goes_imager_projection'].perspective_point_height

# Longitude of the satellite's orbit
lon = g16nc.variables['goes_imager_projection'].longitude_of_projection_origin

# Which direction do we sweep in?
sweep = g16nc.variables['goes_imager_projection'].sweep_angle_axis

X = g16nc.variables['x'][:] * height
Y = g16nc.variables['y'][:] * height

g16nc.close()
g16nc = None

#Lets use Basemap to map state & country borders now.
fig = plt.figure(figsize=(6,6),dpi=200)
m = Basemap(resolution='l', projection='geos', lon_0=lon,
            llcrnrx=X.min(), llcrnry=Y.min(), urcrnrx=X.max(), urcrnry=Y.max())
m.imshow(np.flipud(truecolor))
m.drawcountries()
m.drawstates()
plt.show()

#-----------------------------------------------------------------------------------

#GLM Point Data

#The Geostationary Lightning Mapper (GLM) provides point based data of lightning, flashes,
#groups, and events.
g16glm = Dataset('OR_GLM-L2-LCFA_G16_s20180471253200_e20180471253400_c20180471253551.nc', 'r')

#GLM Counts
#GLM provides latitude longitude points for events (detections) which are then aggregated
#into group and flash events. Think about how for a single ground strike in a supercell
#thunderstorm the intracloud portion of the lightning may extend for many 10s to 100s of
#miles as the extent of the events. This data is then aggregated to provide cluster centroids
#which are stored as the flashes.

event_lat = g16glm.variables['event_lat'][:]
event_lon = g16glm.variables['event_lon'][:]

group_lat = g16glm.variables['group_lat'][:]
group_lon = g16glm.variables['group_lon'][:]

flash_lat = g16glm.variables['flash_lat'][:]
flash_lon = g16glm.variables['flash_lon'][:]

#GLM Plots

#Here we will just plot all of the data on top of each other so you can observed how the data
#are aggregated together.

fig = plt.figure(figsize=(6,6),dpi=200)
map = Basemap(projection='merc', lat_0 = 0, lon_0 = -70.0,
    resolution = 'l', area_thresh = 1000.0,
    llcrnrlon=-135.0, llcrnrlat=-65.0,
    urcrnrlon=0.0, urcrnrlat=65.0)
 
map.drawcoastlines()
map.drawcountries()
map.fillcontinents(color = 'tan')
map.drawmapboundary()


# Plot events as large blue dots
event_x, event_y = map(event_lon, event_lat)
map.plot(event_x, event_y, 'bo', markersize=7)

# Plot groups as medium green dots
group_x, group_y = map(group_lon, group_lat)
map.plot(group_x, group_y, 'go', markersize=3)

# Plot flashes as small red dots
flash_x, flash_y = map(flash_lon, flash_lat)
map.plot(flash_x, flash_y, 'ro', markersize=1)
 
plt.show()

