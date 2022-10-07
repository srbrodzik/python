import numpy as np
import pyart
import os
import glob
import pdb

nc_file = '/home/disk/shear2/brodzik/python/raintype/polar/input/20111016/cfrad.20111016_120032.049_to_20111016_120528.178_SPOL_v4117_SUR.nc'
radar = pyart.io.read(nc_file)
for ray in radar.sweep_number['data']
azimuth = radar.fixed_angle['data'][ray]
