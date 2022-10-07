#!/usr/bin/python3

'''
**Function to compute the distance in km from a distance in degrees for a region of
**the Earth centered at the lat-lon coordinates
deg    - distance measured in degrees
kms[0] - Zonal distance (Equivalent of 1°) in kilometers at given Lon, Lat
kms[1] - meridional distance (Equivalent of 1° )in Kiometers at given Lon, Lat and
Re=6378.137   this is the earth's radius at the equator   
'''
import math
import numpy as np

def deg2km(deg,lon,lat):
    # mean earth radius (this is equivalent to the matlab function)
    Re=6371.
    
    kms=np.zeros(2)
    kms[0]=Re/(180.0/math.pi)
    kms[1]=(Re/(180.0/math.pi))*math.cos(lat*(math.pi/180))
    kms=deg*kms
    return kms

'''
Function to read in topography data and set all alts less than zero to zero
INPUT - topo_file
OUTPUT - topo_lat
         topo.lon
         DEM
'''
def read_topo(topo_file):

    # ONLY WORKS WITH PYTHON2, NOT PYTHON3
    
    import hdf
    import numpy as np
    from pyhdf import SD
    
    # for testing
    topo_file = '/home/disk/shear2/brodzik/GTOPO30/gtopo5km_WMP_gpm_v06.hdf'
    
    lat_name = 'Latitude'
    lon_name = 'Longitude'
    elev_name = 'elevation'

    # open file
    hdf = SD.SD(topo_file)
    
    sds = hdf.select(lat_name)
    topo_lat =sds.get()
    sds.endaccess()
    
    sds = hdf.select(lon_name)
    topo_lon = sds.get()
    sds.endaccess()
    
    sds = hdf.select(elev_name)
    topo_elev = sds.get()
    sds.endaccess()
    topo_elev[topo_elev <= 0] = 0

    num_lats = topo_lat.shape
    num_lons = topo_lon.shape
    # topo_elev.shape = num_lats,num_lons

    hdf.end()
    
    return topo_lat,topo_lon,topo_elev

def read_topo_h5(topo_file):

    import numpy as np
    import h5py
    
    # for testing
    #topo_file = '/home/disk/shear2/brodzik/GTOPO30/gtopo5km_AFC_gpm.h5'
    topo_file = '/home/disk/blitz/data/GTOPO30/gtopo5km_AFC_gpm.h5'
    
    hf = h5py.File(topo_file,'r')
    hf.keys()
    topo_lat = np.array(hf.get('lats005'))
    topo_lon = np.array(hf.get('lons005'))
    topo_elev = np.array(hf.get('topo005'))

    hf.close()

    return topo_lat,topo_lon,topo_elev
    
