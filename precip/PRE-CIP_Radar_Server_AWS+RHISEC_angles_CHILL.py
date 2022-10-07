
# coding: utf-8

# # Using Python to Access NCEI Archived NEXRAD Level 2 Data & Plotting Azimuths w.r.t. Radar Location -- CHILL

# This notebook shows how to access the THREDDS Data Server (TDS) instance that is serving up archived NEXRAD Level 2 data hosted on Amazon S3. The TDS provides a mechanism to query for available data files, as well as provides access to the data as native volume files, through OPeNDAP, and using its own CDMRemote protocol. Since we're using Python, we can take advantage of Unidata's Siphon package, which provides an easy API for talking to THREDDS servers.
# 
# **NOTE:** Due to data charges, the TDS instance in AWS only allows access to .edu domains. For other users interested in using Siphon to access radar data, you can access recent (2 weeks') data by changing the server URL below to: http://thredds.ucar.edu/thredds/radarServer/nexrad/level2/IDD/
# 
# **But first!**
# Bookmark these resources for when you want to use Siphon later!
# + [latest Siphon documentation](http://siphon.readthedocs.org/en/latest/)
# + [Siphon github repo](https://github.com/Unidata/siphon)
# + [TDS documentation](http://www.unidata.ucar.edu/software/thredds/current/tds/TDS.html)

# ## Downloading the single latest volume
# 

# Just a bit of initial set-up to use inline figures and quiet some warnings.

# In[3]:

import matplotlib
import warnings
warnings.filterwarnings("ignore", category=matplotlib.cbook.MatplotlibDeprecationWarning)
get_ipython().magic(u'matplotlib inline')


# First we'll create an instance of RadarServer to point to the appropriate radar server access URL.

# In[ ]:

from siphon.radarserver import RadarServer
rs = RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')


# Next, we'll create a new query object to help request the data. Using the chaining methods, let's ask for the latest data at the radar KCYS (Cheyenne, WY). We see that when the query is represented as a string, it shows the encoded URL.

# In[ ]:

from datetime import datetime, timedelta
query = rs.query()
query.stations('KCYS').time(datetime.utcnow())


# We can use the RadarServer instance to check our query, to make sure we have required parameters and that we have chosen valid station(s) and variable(s)
# 

# In[ ]:

rs.validate_query(query)


# Make the request, which returns an instance of TDSCatalog; this handles parsing the returned XML information.

# In[ ]:

catalog = rs.get_catalog(query)


# We can look at the datasets on the catalog to see what data we found by the query. We find one volume in the return, since we asked for the volume nearest to a single time.

# In[ ]:

catalog.datasets


# We can pull that dataset out of the dictionary and look at the available access URLs. We see URLs for OPeNDAP, CDMRemote, and HTTPServer (direct download).

# In[ ]:

ds = list(catalog.datasets.values())[0]
ds.access_urls


# We'll use the CDMRemote reader in Siphon and pass it the appropriate access URL.

# In[ ]:

from siphon.cdmr import Dataset
data = Dataset(ds.access_urls['CdmRemote'])


# We define some helper functions to make working with the data easier. One takes the raw data and converts it to floating point values with the missing data points appropriately marked. The other helps with converting the polar coordinates (azimuth and range) to Cartesian (x and y).

# In[ ]:

import numpy as np
def raw_to_masked_float(var, data):
    # Values come back signed. If the _Unsigned attribute is set, we need to convert
    # from the range [-127, 128] to [0, 255].
    if var._Unsigned:
        data = data & 255

    # Mask missing points
    data = np.ma.array(data, mask=data==0)

    # Convert to float using the scale and offset
    return data * var.scale_factor + var.add_offset

def polar_to_cartesian(az, rng):
    az_rad = np.deg2rad(az)[:, None]
    x = rng * np.sin(az_rad)
    y = rng * np.cos(az_rad)
    return x, y


# The CDMRemote reader provides an interface that is almost identical to the usual python NetCDF interface. We pull out the variables we need for azimuth and range, as well as the data itself.

# In[ ]:

sweep = 0
ref_var = data.variables['Reflectivity_HI']
ref_data = ref_var[sweep]
rng = data.variables['distanceR_HI'][:]
az = data.variables['azimuthR_HI'][sweep]


# Then convert the raw data to floating point values and the polar coordinates to Cartesian.

# In[ ]:

ref = raw_to_masked_float(ref_var, ref_data)
x, y = polar_to_cartesian(az, rng)


# MetPy is a Python package for meteorology (Documentation: http://metpy.readthedocs.org and GitHub: http://github.com/MetPy/MetPy). We import MetPy and use it to get the colortable and value mapping information for the NWS Reflectivity data.

# In[ ]:

from metpy.plots import ctables  # For NWS colortable
ref_norm, ref_cmap = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20, 0.5)


# Finally, we plot them up using matplotlib and cartopy. We create a helper function for making a map to keep things simpler later.

# In[ ]:

import matplotlib.pyplot as plt
import cartopy
import cartopy.io.shapereader as shpreader
import cartopy.geodesic as cgds
#increase font size
plt.rcParams.update({'font.size': 18})
import matplotlib.patches as mpatches
import shapely
from pyproj import Geod

def new_map(fig, lon, lat, azimuths):
    # Create projection centered on the radar. This allows us to use x
    # and y relative to the radar.
    proj = cartopy.crs.LambertConformal(central_longitude=lon, central_latitude=lat)

    # New axes with the specified projection
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    # Add coastlines
    ax.coastlines('10m', 'black', linewidth=2, zorder=2)
    
    #RRB: zoom in (or at least try to...)
    ax.set_extent([lon-2, lon+2, lat-2.5, lat+1])
    
    # Grab state borders
    state_borders = cartopy.feature.NaturalEarthFeature(
        category='cultural', name='admin_1_states_provinces_lines',
        scale='50m', facecolor='none')
    ax.add_feature(state_borders, edgecolor='black', linewidth=1, zorder=3)
    
    #add county borders
    reader = shpreader.Reader('countyl010g.shp')
    counties = list(reader.geometries())
    COUNTIES = cartopy.feature.ShapelyFeature(counties, cartopy.crs.PlateCarree())
    ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')
    
    #add star at the location of the CHIVO radar
    lat_CHIVO=40.44625
    lon_CHIVO=-104.63708
    ax.scatter(lon_CHIVO,lat_CHIVO,marker='*',s=100,color='k',transform=cartopy.crs.PlateCarree())

    #add a ring at ~50 km from CHIVO
    circle_points = cgds.Geodesic().circle(lon=lon_CHIVO, lat=lat_CHIVO, radius=50000, n_samples=100, endpoint=False)
    geom = shapely.geometry.Polygon(circle_points)
    ax.add_geometries((geom,), crs=cartopy.crs.PlateCarree(), facecolor='none', edgecolor='red', linewidth=2)

    #add a ring at ~75 km from CHIVO
    #circle_points = cgds.Geodesic().circle(lon=lon_CHIVO, lat=lat_CHIVO, radius=75000, n_samples=100, endpoint=False)
    #geom = shapely.geometry.Polygon(circle_points)
    #ax.add_geometries((geom,), crs=cartopy.crs.PlateCarree(), facecolor='none', edgecolor='red', linewidth=2)
    
    #add a ring at ~100 km from CHIVO
    circle_points = cgds.Geodesic().circle(lon=lon_CHIVO, lat=lat_CHIVO, radius=100000, n_samples=100, endpoint=False)
    geom = shapely.geometry.Polygon(circle_points)
    ax.add_geometries((geom,), crs=cartopy.crs.PlateCarree(), facecolor='none', edgecolor='red', linewidth=2)

    #add a ring at ~150 km from CHIVO
    circle_points = cgds.Geodesic().circle(lon=lon_CHIVO, lat=lat_CHIVO, radius=150000, n_samples=100, endpoint=False)
    geom = shapely.geometry.Polygon(circle_points)
    ax.add_geometries((geom,), crs=cartopy.crs.PlateCarree(), facecolor='none', edgecolor='red', linewidth=2)
    
    #line points at azimuths
    geod = Geod(ellps='sphere')
    rng = 150 * 1000
    for az in azimuths:
        lon2, lat2, _ = geod.fwd(lon_CHIVO, lat_CHIVO, az, rng)
        ax.plot([lon_CHIVO,lon2],[lat_CHIVO,lat2],'--k',transform=cartopy.crs.PlateCarree())
        ax.text(lon2,lat2,az,transform=cartopy.crs.PlateCarree())
    
    return ax


# ##### Use the function to make a new map and plot a colormapped view of the data

# In[ ]:

fig = plt.figure(figsize=(20, 20))

az = np.arange(137, 200, 15)
#az = []
#az = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350]
ax = new_map(fig, data.StationLongitude, data.StationLatitude, az)
cf = ax.pcolormesh(x, y, ref, cmap=ref_cmap, norm=ref_norm, zorder=0)
cbar_ax = fig.add_axes([0.92, 0.225, 0.015, 0.55])
cbar = fig.colorbar(cf, cax=cbar_ax)
plt.title(datetime.utcnow())
plt.savefig('KYCS+CHILL_with100kmring_and_azimuths.png',bbox_inches='tight')


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



