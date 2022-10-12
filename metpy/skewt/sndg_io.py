import os
import sys
import numpy as np
import pandas as pd
import xarray as xr
import datetime
from calculations import calculate_dewpoint

MPS2KTS = 1.94
col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']

def read_infile(inpath,infile,format):

    file_out_dt_fmt = '%Y%m%d%H%M'
    title_dt_fmt = '%Y%m%d %H%MZ'

    if format == 'UIUCnc':
        ds = xr.open_dataset(inpath+'/'+infile)

        # determine out_fname and figtitle
        attrs = ds.attrs
        stn_id = 'UILL_Mobile'
        file_time = datetime.datetime.strptime(attrs['start_datetime'],'%Y-%m-%dT%H:%M:%SZ')
        location_str = attrs['location']
        location_list = location_str.split()
        lat = float(location_list[0])
        lon = float(location_list[3])
        if location_list[-1] == 'west':
            lon *= -1
        if location_list[2] == 'south':
            lat *= -1
        out_fname = 'upperair.UILL_sonde.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
        figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)

        # read input data into dataframe object
        temp = ds['TC'].values
        height = ds['HAGL'].values
        # Joe Finlon fix for m/s to kts
        wspd = ds['WINDSPD'].values*MPS2KTS
        pres = ds['PRESS'].values
        drct = ds['WINDDRN'].values
        RH = ds['RH'].values
        dew = []
        for i,val in enumerate(RH):
            dew.append(calculate_dewpoint(temp[i], RH[i]))
        # dewpt needs to be an np.array, not a list, so convert it
        dewpt = np.asarray(dew)

    else:
        print('test')

    # Create dataframe
    df = pd.DataFrame({'pressure':pres,
                       'height':height,
                       'temperature':temp,
                       'dewpoint':dewpt,
                       'direction':drct,
                       'speed':wspd})
    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'),
                   how='all').reset_index(drop=True)

    return df, out_fname, figtitle

