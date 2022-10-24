import os
import re
import io
import sys
import numpy as np
import pandas as pd
import xarray as xr
import datetime
from calculations import calculate_dewpoint

MPS2KTS = 1.94
col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']
file_out_dt_fmt = '%Y%m%d%H%M'
title_dt_fmt = '%Y%m%d %H%MZ'
file_in_dt_fmt_UWYO = '%Y%m%d%H%M'
file_in_dt_fmt_CSU = '%Y%m%d_%H%M'

def read_infile(inpath,infile,fmt):
    if fmt == 'Albany':
        (df,out_fname,figtitle) = read_Albany(inpath,infile)
    elif fmt == 'CSU':
        (df,out_fname,figtitle) = read_CSU(inpath,infile)
    elif fmt == 'MUtxt_ws':
        (df,out_fname,figtitle) = read_MUtxt_ws(inpath,infile)
    elif fmt in ['NCSU','Purdue']:
        (df,out_fname,figtitle) = read_sharppy(inpath,infile,fmt)
    elif fmt == 'RTSO':
        (df,out_fname,figtitle) = read_RTSO(inpath,infile)
    elif fmt in ['SBUnc_mobile','SBUnc','MUnc']:
        (df,out_fname,figtitle) = read_nc_file(inpath,infile,fmt)
    elif fmt == 'UIUCnc':
        (df,out_fname,figtitle) = read_UIUCnc(inpath,infile)
    elif fmt == 'UNCA':
        (df,out_fname,figtitle) = read_UNCA(inpath,infile)
    elif fmt == 'UWYO':
        (df,out_fname,figtitle) = read_UWYO(inpath,infile)
    elif fmt == 'Valpo':
        (df,out_fname,figtitle) = read_Valpo(inpath,infile)
    else:
        print('Unrecognized fmt = ',fmt)
        df = pd.DataFrame()
        out_fname = ''
        figtitle = ''
    return df, out_fname, figtitle

def read_Albany(inpath,infile):

    # parse file name
    (base,ext) = os.path.splitext(infile)
    (category,platform,datetimeStr,product) = base.split('.')
    (stn_id_plot,location) = product.split('-')

    stn_id_plot = 'U'+stn_id_plot
    stn_id_file = 'UALB'

    # get file time
    file_time = datetime.datetime.strptime(datetimeStr,'%Y%m%d%H%M')

    # open file and read in column headers; replace names as required by calling code
    with open(inpath+'/'+infile, encoding="utf8", errors='ignore') as f:
        line = f.readline()
        line = line.strip()
        columns = line.split()
    columns = [sub.replace('Alt_AGL', 'height') for sub in columns]
    columns = [sub.replace('Press', 'pressure') for sub in columns]
    columns = [sub.replace('Temp', 'temperature') for sub in columns]
    columns = [sub.replace('DP', 'dewpoint') for sub in columns]
    columns = [sub.replace('WDirn', 'direction') for sub in columns]
    columns = [sub.replace('WSpeed', 'speed') for sub in columns]

    # read data into dataframe
    df = pd.read_csv(inpath+'/'+infile, skiprows=3, encoding="latin1", delim_whitespace=True, header=None)
    columns.insert(2,'AM-PM')
    df.columns = columns

    # get lat and lon
    latStr = df.iloc[0]['Lat/N']
    deg = int(latStr[0:2])
    min = int(latStr[3:5])
    sec = float(latStr[6:10])
    dir = latStr[-1]
    lat = dirMinSec2deg(deg,min,sec,dir)

    lonStr = df.iloc[0]['Long/E']
    deg = int(lonStr[0:3])
    min = int(lonStr[4:6])
    sec = float(lonStr[7:11])
    dir = lonStr[-1]
    lon = dirMinSec2deg(deg,min,sec,dir)

    # define out_fname and figtitle
    out_fname = 'upperair.{id}_sonde.{dt}.skewT.png'.format(id = stn_id_file, dt = file_time.strftime(file_out_dt_fmt))
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id_plot, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)  
    
    for key in df.keys():
	# if missing values in data, have to replace with NaNs and then convet to numeric as it
        # originaly is read strings instead of float64; coerce flag nicely takes anything that
        # cannot be converted and puts NaN
        if df[key].dtype.name != 'float64':
            df[key] = pd.to_numeric(df[key], errors='coerce')

    # convert wind speed from m/s to knots
    df['speed'] = df['speed'] * MPS2KTS

    return df, out_fname, figtitle

def read_CSU(inpath,infile):

    rawfile = open(inpath+'/'+infile, encoding="latin-1").readlines()

    # I don't think we need this try/except loop
    try:
        serial_number = rawfile[3].split('\t')[1].strip()
        station = '%s, %s'%(rawfile[0].split('\t')[1].strip(), serial_number)
        time_string = ' '.join(rawfile[5].split('\t')[1:]).strip()
        sounding_date = datetime.datetime.strptime(time_string, '%Y-%m-%dT%H:%M:%S')

    except Exception:
        pass
        station = 'test'
        sounding_date = None

    # Find the Header line by looping thru the file lines
    for il, line in enumerate(rawfile):
        if 'elapsed' in line.lower():
            header_line = il
            break
        else:
            pass

    header_split = rawfile[header_line].lower().split()

    htcol = header_split.index('heightmsl')
    prescol = header_split.index('p')
    tempcol = header_split.index('temp')
    dewcol = header_split.index('dewp')
    spdcol = header_split.index('speed')
    drctcol = header_split.index('dir')
    latcol = header_split.index('lat')
    loncol = header_split.index('lon')
    
    if 'time' in header_split:
        htcol -= 1; prescol -= 1; tempcol -= 1; dewcol -= 1; spdcol -= 1; drctcol -= 1; latcol -= 1; loncol -= 1

    file_data = np.genfromtxt(inpath+'/'+infile, encoding="latin-1", skip_header=header_line+2)

    height = file_data[:,htcol]
    pres = file_data[:,prescol]
    temp = file_data[:,tempcol]
    dew = file_data[:,dewcol]
    # Convert wspd from m/s to knots
    wspd = file_data[:,spdcol]*MPS2KTS
    drct = file_data[:,drctcol]
    lon = file_data[:,loncol]
    lat = file_data[:,latcol]

    drct[drct > 360] = 0
    wspd[wspd > 999] = 0

    df = pd.DataFrame({'pressure':pres,
                       'height':height,
                       'temperature':temp,
                       'dewpoint':dew,
                       'direction':drct,
                       'speed':wspd})

    # Get out_fname and figtitle
    lat_grnd = lat[0]
    lon_grnd = lon[0]
    file_time_string = infile[4:17]
    file_time = datetime.datetime.strptime(file_time_string, file_in_dt_fmt_CSU)
    out_fname = 'upperair.CSU_sonde.{dt}.skewT.png'.format(dt=file_time.strftime(file_out_dt_fmt))
    figtitle = 'CSU {dt} sounding ({lat:.3f} {lon:.3f})'.format(dt=file_time.strftime(title_dt_fmt), lat=lat_grnd, lon=lon_grnd)

    return df, out_fname, figtitle

def read_MUtxt_ws(inpath,infile):

    #### NEED TO ADD LAT AND LON TO FIGURE TITLE ####
    
    # determine out_fname and figtitle
    # Assumes infile naming convention is: upperair.UMILL_windsonde1.202202250600.txt
    parts = infile.split('.')
    snd_id = parts[1]
    file_time_string = parts[2]
        
    file_time = datetime.datetime.strptime(file_time_string, file_in_dt_fmt_UWYO)
    out_fname = 'upperair.{snd}.{dt}.skewT.png'.format(snd=snd_id, dt=file_time.strftime(file_out_dt_fmt))
    figtitle = '{snd} {dt} sounding'.format(snd=snd_id, dt=file_time.strftime(title_dt_fmt))

    # read in file after header as DataFrame
    df = pd.read_csv(inpath+'/'+infile, skiprows=7, header=None)

    # Create column headings
    df.columns = ['pressure','height','temp','dewpt','wdir','wspd']

    # Remove last row ('%END%')
    while df.iloc[-1]['pressure'] == '%END%':
        numRows = len(df.index)
        df.drop([numRows-1],inplace=True)

    for key in df.keys():
	# if missing values in data, have to replace with NaNs and then convet to numeric as it
        # originaly is read strings instead of float64; coerce flag nicely takes anything that
        # cannot be converted and puts NaN
        if df[key].dtype.name != 'float64':
            df[key] = pd.to_numeric(df[key], errors='coerce')
        if 'pressure' in key:
            pres = (df[key][0:])
        elif 'height' in key:
            height = (df[key][0:])
        elif 'temp' in key:
            temp = (df[key][0:])
        elif 'dewpt' in key:
            dewpt = (df[key][0:])
        elif 'wdir' in key:
            drct = (df[key][0:])
        elif 'wspd' in key:
            # this format already in knots; no need to convert
            wspd = (df[key][0:])

    # Create dataframe
    df = pd.DataFrame({'pressure':pres,
                       'height':height,
                       'temperature':temp,
                       'dewpoint':dewpt,
                       'direction':drct,
                       'speed':wspd})

    #### NEED TO CHANGE -9999.00 VALUES TO NaN ####
    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'),
                   how='all').reset_index(drop=True)
        
    return df, out_fname, figtitle

def read_sharppy(inpath,infile,fmt):

    # Get out_fname and figtitle
    if fmt == 'NCSU':
        stn_id = 'NCSU'
    elif fmt == 'Purdue':
        stn_id = 'PurdueU'

    # read file header
    with open(inpath+'/'+infile) as myfile:
        head = myfile.readline()
    parts = head.split(' ')
    for part in parts:
        if 'lat' in part:
            (junk,latStr) = part.split('=')
        elif 'lon' in part:
            (junk,lonStr) = part.split('=')
        elif 'time' in part:
            (junk,dateStr) = part.split('=')
    timeStr = parts[-1].rstrip()
                
    # get file time
    datetimeStr = (dateStr+'T'+timeStr)
    file_time_est = datetime.datetime.strptime(datetimeStr,'%Y-%m-%dT%H:%M')
    file_time_utc = file_time_est + datetime.timedelta(hours=5)

    # turn lat/lon into floats
    lat = float(latStr[:-1])
    lon = float(lonStr[:-1])

    # define out_fname and figtitle
    out_fname = 'upperair.{fmt}_sonde.{dt}.skewT.png'.format(fmt = fmt.upper(), dt = file_time_utc.strftime(file_out_dt_fmt))
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id, dt = file_time_utc.strftime(title_dt_fmt),
                                                                     lati=lat, long=lon)
        
    # Create dataframe and rename columns
    df = pd.read_csv(inpath+'/'+infile, skiprows=1)
    df.columns = ['height','pressure','temperature','rh','speed','direction']

    # Convert speed from m/s to knots
    df['speed'] = df['speed'] * MPS2KTS

    # Handle missing values in data
    for key in df.keys():
	# Replace missing values with NaNs and then convert to numeric as they are originally read 
        # as strings instead of float64; coerce flag nicely takes anything that cannot be converted
        # and substitutes NaN
        if df[key].dtype.name != 'float64':
            df[key] = pd.to_numeric(df[key], errors='coerce')

    # Calculate dewpoint and add to dataframe
    dew = []
    for i,val in enumerate(df['rh']):
        dew.append(calculate_dewpoint(df['temperature'][i], df['rh'][i]))
    dewpt = np.asarray(dew)
    df.insert(3, "dewpoint", dewpt, True)
    
    return df, out_fname, figtitle

def read_RTSO(inpath,infile):

    # read file header
    with open(inpath+'/'+infile) as myfile:
        #head = [next(myfile) for x in xrange(8)]  # for python2
        head = [next(myfile) for x in range(8)]  # for python3: xrange changed to range
    dateStr = head[2].strip()
    stationStr = head[3].strip()
    latlonStr = head[4].strip()

    # get file time
    temp = dateStr.split(' ')
    datetimeStr = temp[1]+' '+temp[4]
    file_time = datetime.datetime.strptime(datetimeStr,'%m/%d/%Y %H:%M:%S')

    # get station id
    stn_id = stationStr
    if stn_id == 'N-179':
        stn_id_file = 'Wallops_VA'
        stn_id_plot = 'WAL'
    else:
        stn_id_file = 'unknown'
        stn_id_plot = 'UNK'
        
    # get lat/lon
    temp = latlonStr.split(' ')
    lat = abs(float(temp[2][:-1]))
    if temp[2][-1] == 'S':
        lat = -lat
    lon = abs(float(temp[4][:-1]))
    if temp[4][-1] == 'W':
        lon = -lon

    # define out_fname and figtitle
    out_fname = 'upperair.SkewT.{dt}.{stn}.png'.format(dt = file_time.strftime(file_out_dt_fmt), stn = stn_id_file)
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id_plot, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)

    # read in file after header as DataFrame
    df = pd.read_csv(inpath+'/'+infile, skiprows=9, delim_whitespace=True)

    #df.columns = ['time (s)','pres (mb)','temp (degC)','rel_hum (%)','geop_ht (m)','dew_pt (degC)',
    #              'ref_ind (N units)','grad_ref_ind (N/km)','mod_ref_ind (M units)','speed_snd (m/s)',
    #              'air_dens (g/m3)','vap_pres (mb)','pot_temp (degC)','vir_temp (degC)','spec_hum (g/kg)',
    #              'spare1','spare2','spare3','spare4','utc_time (s)','wspd (m/s)','wdir (deg)',
    #              'ns_wind_comp (m/s)','ew_wind_comp (m/s)','vert_wind_comp (m/s)','lon (deg)','lat (deg)',
    #              'geom_ht (m)']

    df.columns = ['time (s)','pressure','temperature','rel_hum (%)','geop_ht (m)','dewpoint',
                  'ref_ind (N units)','grad_ref_ind (N/km)','mod_ref_ind (M units)','speed_snd (m/s)',
                  'air_dens (g/m3)','vap_pres (mb)','pot_temp (degC)','vir_temp (degC)','spec_hum (g/kg)',
                  'spare1','spare2','spare3','spare4','utc_time (s)','speed','direction',
                  'ns_wind_comp (m/s)','ew_wind_comp (m/s)','vert_wind_comp (m/s)','lon (deg)','lat (deg)',
                  'height']

    # Handle missing values in data
    for key in df.keys():
	# Replace missing values with NaNs and then convert to numeric as they are originally read 
        # as strings instead of float64; coerce flag nicely takes anything that cannot be converted
        # and substitutes NaN
        if df[key].dtype.name != 'float64':
            df[key] = pd.to_numeric(df[key], errors='coerce')

    # Remove rows with pressure val = NaN; we need this for plotting wind barbs
    df = df[df['pressure'].notna()]
    df = df[df['dewpoint'].notna()]
            
    # Convert wind speed from m/s to knots
    df['speed'] = df['speed'] * MPS2KTS
    
    return df, out_fname, figtitle

def read_nc_file(inpath,infile,fmt):

    # The data for SBUnc_mobile, SBUnc and MUnc is in the same format so
    # we use this to read it and fill the dataframe

    ds = xr.open_dataset(inpath+'/'+infile,decode_times=False)

    # get attributes needed for out_fname and figtitle
    attrs = ds.attrs
    file_time = datetime.datetime.strptime(attrs['ReleaseTime'],'%Y/%m/%d %H:%M:%S')
    lat_str_raw = attrs['SiteLocation_Latitude']
    lon_str_raw = attrs['SiteLocation_Longitude']
    (lat_str,lat_dir) = lat_str_raw.split()
    (lon_str,lon_dir) = lon_str_raw.split()
    lat = float(lat_str)
    lon = float(lon_str)

    if fmt == 'SBUnc_mobile':
        stn_id = 'SBU_Mobile'
        out_fname = 'upperair.SBU_sonde.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    elif fmt == 'SBUnc':
        stn_id = 'SBU'
        out_fname = 'upperair.SBU_sonde.{dt}.Stonybrook_NY_skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    elif fmt == 'MUnc':
        stn_id = 'UMILL'
        out_fname = 'upperair.UMILL_sonde.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)
            
    # read input data into dataframe object
    # Data tends to oscillate, so end index is the first time the heights begin to decrease
    # All data is taken if heights never decrease
    # try:
    #    end_index = np.where(np.gradient(ds['geometric_height'].values) < 0  )[0][0]
    # except:
    #    end_index = -1
    end_index = -1
    pres = ds['pressure'].values[:end_index]
    height = ds['geometric_height'].values[:end_index]
    temp = ds['temperature'].values[:end_index]
    dewpt = ds['dewpoint_temperature'].values[:end_index]
    drct = ds['wind_direction'].values[:end_index]
    wspd = ds['wind_speed'].values[:end_index]*MPS2KTS

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

def read_UIUCnc(inpath,infile):

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
    pres = ds['PRESS'].values
    height = ds['HAGL'].values
    temp = ds['TC'].values
    RH = ds['RH'].values
    dew = []
    for i,val in enumerate(RH):
        dew.append(calculate_dewpoint(temp[i], RH[i]))
    # dewpt needs to be an np.array, not a list, so convert it
    dewpt = np.asarray(dew)
    drct = ds['WINDDRN'].values
    wspd = ds['WINDSPD'].values*MPS2KTS

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

def read_UNCA(inpath,infile):

    # read file header
    #with open(inpath+'/'+infile, encoding="latin1") as myfile:
    with io.open(inpath+'/'+infile, encoding="cp1252") as myfile:
        #head = [next(myfile) for x in xrange(8)]  # for python2
        head = [next(myfile) for x in range(22)]  # for python3: xrange changed to range
    dateStr = head[-1].strip()
    stationStr = head[3].strip()
    latStr = head[7].strip()
    lonStr = head[8].strip()
    #columns = head[50].strip().split()

    # get file time
    temp = dateStr.split()
    datetimeStr = temp[3]+' '+temp[4]
    AmPm = temp[5]
    file_time = datetime.datetime.strptime(datetimeStr,'%m/%d/%Y %H:%M:%S')
    if AmPm == 'PM':
        file_time = file_time + datetime.timedelta(hours=12)

    # get station id
    temp = stationStr.split(':')
    stn_id_file = temp[1].strip().replace(' ','_')
    stn_id_plot = 'UNCA'
        
    # get lat/lon
    temp = latStr.split()
    (deg,min,sec,dir) = re.sub(r"[^a-zA-Z0-9\. ]", " ", temp[-1]).split()
    lat = dirMinSec2deg(float(deg),float(min),float(sec),dir)

    temp = lonStr.split()
    (deg,min,sec,dir) = re.sub(r"[^a-zA-Z0-9\. ]", " ", temp[-1]).split()
    lon = dirMinSec2deg(float(deg),float(min),float(sec),dir)

    # define out_fname and figtitle
    out_fname = 'upperair.{stn}_sonde.{dt}.skewT.png'.format(stn = stn_id_file, dt = file_time.strftime(file_out_dt_fmt))
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id_file, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)

    # find data for plotting past header and significant levels
    with io.open(inpath+'/'+infile, encoding="cp1252") as f:
        for num,line in enumerate(f,1):
            newline = line.strip()
            if newline.startswith('FltTime'):
                linenum = num
                #print(linenum)
                break

    with io.open(inpath+'/'+infile, encoding="cp1252") as f:
        lines = f.readlines()
        data = lines[linenum+2:]
        fout = open('/tmp/sndg.dat','w')
        for num,line in enumerate(data):
            fout.write(data[num])
        fout.close()
        
    df = pd.read_csv('/tmp/sndg.dat', skiprows=0, encoding='latin1', delim_whitespace=True, header=None)
    df.columns = (['FltTime','pressure','temperature','rh','speed','direction','height'])

    # add dewpoint to dataframe
    dew = []
    for i,val in enumerate(df['rh']):
        dew.append(calculate_dewpoint(df['temperature'][i], df['rh'][i]))
    dewpt = np.asarray(dew)
    df.insert(3, "dewpoint", dewpt, True)

    # convert wind speed from m/s to knots
    df['speed'] = df['speed'] * MPS2KTS

    os.remove('/tmp/sndg.dat')
        
    return df, out_fname, figtitle

def read_UWYO(inpath,infile):

    # Get out_fname and figtitle
    parts = infile.split('.')
    file_time_string = parts[2]
    snd_id = parts[3]
    
    file_time = datetime.datetime.strptime(file_time_string, file_in_dt_fmt_UWYO)
    out_fname = 'upperair.SkewT.{dt}.{snd}.png'.format(dt=file_time.strftime(file_out_dt_fmt),snd=snd_id)
    figtitle = '{snd} {dt} sounding'.format(snd=snd_id, dt=file_time.strftime(title_dt_fmt))
    
    # Read data into dataframe
    fid = open(inpath+'/'+infile)
    lines = [line for line in fid.readlines() if line.strip()]
    nlines = len(lines)
    ndata = nlines-34
    output = {}

    fields = lines[3].split()
    units = lines[4].split()

    # First line for WRF profiles differs from the UWYO soundings
    header = lines[1]
    if header[:5] == '00000':
        # WRF profile
        station = '-99999'
        Longitude = float(header.split()[5].strip(","))
        Latitude = float(header.split()[6])
        sounding_date = header.split()[-1]
    else:
        station = header[:5]
        dstr = (' ').join(header.split()[-4:])
        sounding_date = datetime.datetime.strptime(dstr, "%HZ %d %b %Y").strftime("%Y-%m-%d_%H:%M:%S") 

        #if station_name is not None: station = station_name
            
        file_data = np.genfromtxt(inpath+'/'+infile, skip_header=18, skip_footer=51)
           
        htcol = fields.index('HGHT')
        prescol = fields.index('PRES')
        tempcol = fields.index('TEMP')
        dewcol = fields.index('DWPT')
        spdcol = fields.index('SKNT')
        drctcol = fields.index('DRCT')
            
        height = file_data[:,htcol]
        pres = file_data[:,prescol]
        temp = file_data[:,tempcol]
        dewpt = file_data[:,dewcol]
        # wspd data in kts
        wspd = file_data[:,spdcol]
        drct = file_data[:,drctcol]

        drct[drct > 360] = 0
        wspd[wspd > 999] = 0

    # Close input file
    fid.close()
        
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

def read_Valpo(inpath,infile):

    # Get out_fname and figtitle
    latFound = False
    lonFound = False
    with open(inpath+'/'+infile, 'r') as fp:
        while not latFound or not lonFound:
            line = fp.readline().strip()
            #print(line)
            if 'lat' in line:
                lat = float(line.replace('lat: ',''))
                latFound = True
            elif 'lon' in line:
                lon = float(line.replace('lon: ',''))
                lonFound = True

    # Assumes parts of input file name separated by periods & filename = upperair.VALPO_sonde.<datetime>.csv
    parts = infile.split('.')
    snd_id = parts[1]
    file_time_string = parts[2]
        
    file_time = datetime.datetime.strptime(file_time_string, file_in_dt_fmt_UWYO)
    out_fname = 'upperair.{snd}.{dt}.skewT.png'.format(snd=snd_id, dt=file_time.strftime(file_out_dt_fmt))
    figtitle = '{snd} {dt} sounding ({lati:.3f}, {long:.3f})'.format(snd=snd_id, dt=file_time.strftime(title_dt_fmt), lati=lat, long=lon)

    # Create dataframe
    # read in file after header as DataFrame and sort into descending pressures
    # need to do some playing around due to using two sets of significant levels
    df = pd.read_csv(inpath+'/'+infile, skiprows=4, delimiter='\s+')
    df.sort_values(by='Press',ascending=False,inplace=True)
    df = df.drop(labels='Type',axis=1)
    df.reset_index(drop=True,inplace=True)
            
    for key in df.keys():
	# if missing values in data, have to replace with NaNs and then convet to numeric as it
        # originaly is read strings instead of float64; coerce flag nicely takes anything that
        # cannot be converted and puts NaN
        if df[key].dtype.name != 'float64':
            df[key] = pd.to_numeric(df[key], errors='coerce')
        if 'GPM_AGL' in key:
            height = (df[key][0:])
        elif 'Press' in key:
            pres = (df[key][0:])
        elif 'Temp' in key:
            temp = (df[key][0:])
        elif 'DP' in key:
            dewpt = (df[key][0:])
        elif 'WSpeed' in key:
            # wspd in knots
            wspd = (df[key][0:])
        elif 'WDirn' in key:
            drct = (df[key][0:])

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

def dirMinSec2deg(deg,min,sec,dir):
    val = deg+((min+(sec/60))/60)
    if dir == 'W' or dir == 'S':
        val = -val
    return val



