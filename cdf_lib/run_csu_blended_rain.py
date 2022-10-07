"""
run_csu_blended_rain.py

# Stacy Brodzik, UW, May 2016
# brodzik@uw.edu

# wrapper for csu_blended_rain_gan.py for using CSU blended algorithm (see
# Cifelli (2011) for algorithm details)

"""

from __future__ import absolute_import
from __future__ import division
from csu_radartools import csu_blended_rain as cbr
import os
import numpy as np
import warnings
import netCDF4 as nc4

# ---------------------------INPUTS--------------------------
reflDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf/'
rtDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_type_sur/'
outDir = 
