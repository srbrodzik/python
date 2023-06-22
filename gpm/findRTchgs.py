import os
import netCDF4 as nc4
import numpy as np
import datetime as dt
from matplotlib.dates import DateFormatter, MonthLocator

STRA = 1
CONV = 2

#User inputs
indir = '/home/storm/brodzik/Data/CSH'

for fname in os.listdir(indir):
    ncid = nc4.Dataset(indir+'/'+fname,'r')
    rt = np.array(ncid.variables['rain_type'])
    rt_uw = np.array(ncid.variables['rain_type_uw'])
    ncid.close()

    c = np.logical_and(rt==CONV,rt_uw==STRA)
    out = np.column_stack(np.where(c))

