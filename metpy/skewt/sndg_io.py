import os
import re
import io
import sys
import numpy as np
import pandas as pd
import xarray as xr
import datetime
sys.path.append('/home/disk/meso-home/brodzik/python/metpy/skewt')
from calculations import calculate_dewpoint

# NOTE: Oswego data comes in SHARPpy format but without lat/lon in header (use read_sharppy2 to read)
#       If you add lat/lon info to header and format like MUtxt_ws data, can use that input routine
#       Can also use read_sharppy and hard code lat and lon as for NCSU and Purdue
# NOTE: read_sharppy() has launch lat/lon hardwired in the code.  Need to make these optional inputs.

MPS2KTS = 1.94
KM2M = 1000
col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']
file_out_dt_fmt = '%Y%m%d%H%M'
title_dt_fmt = '%Y%m%d %H%MZ'
file_in_dt_fmt_UWYO = '%Y%m%d%H%M'
file_in_dt_fmt_CSU = '%Y%m%d_%H%M'

siteList = {'ABR':'Aberdeen_SD',
            'ALB':'Albany_NY',
            'APX':'Gaylord_MI',
            'BUF':'Buffalo_NY',
            'CAR':'Caribou_ME',
            'CHH':'Chatham_MA',
            'CHS':'Charleston_SC',
            'DTX':'Detroit_MI',
            'DVN':'Davenport_IA',
            'FFC':'Peachtree_City_GA',
            'GRB':'Green_Bay_WI',
            'GSO':'Greensboro_NC',
            'GYX':'Gray_ME',
            'IAD':'Sterling_VA',
            'ILN':'Wilmington_OH',
            'ILX':'Lincoln_IL',
            'INL':'International_Falls_MN',
            'MHX':'Newport_NC',
            'MPX':'Minneapolis_MN',
            'OKX':'Upton_NY',
            'PIT':'Pittsburgh_PA',
            'RNK':'Blacksburg_VA',
            'TBW':'Tampa_Bay_FL',
            'WAL':'Wallops_VA'}

def read_infile(inpath,infile,fmt,lat,lon):
    if fmt in ['Albany_mobile','Albany']:
        (df,out_fname,figtitle) = read_Albany(inpath,infile,fmt)
    elif fmt in ['CMICH','UNCA']:
        (df,out_fname,figtitle) = read_TSPOTINT(inpath,infile,fmt,lat,lon)
    elif fmt == 'CSU':
        (df,out_fname,figtitle) = read_CSU(inpath,infile)
    elif fmt in ['MUtxt_ws','Oswego']:
        (df,out_fname,figtitle) = read_MUtxt_ws(inpath,infile)
    elif fmt in ['NCSU','Purdue','UNDws','HsinChu']:
        (df,out_fname,figtitle) = read_sharppy(inpath,infile,fmt,lat,lon)
    #elif fmt in ['Oswego']:
    #    (df,out_fname,figtitle) = read_sharppy2(inpath,infile)
    elif fmt == 'RTSO':
        (df,out_fname,figtitle) = read_RTSO(inpath,infile)
    elif fmt in ['SBUnc_mobile','SBUnc','MUnc']:
        (df,out_fname,figtitle) = read_nc_file(inpath,infile,fmt)
    elif fmt == 'UIUCnc':
        (df,out_fname,figtitle) = read_UIUCnc(inpath,infile)
    elif fmt == 'UND':
        (df,out_fname,figtitle) = read_gfk(inpath,infile)
        #(df,out_fname,figtitle) = read_gfk_csv(inpath,infile)
    elif fmt == 'UWYO':
        (df,out_fname,figtitle) = read_UWYO(inpath,infile)
    elif fmt == 'Valpo':
        (df,out_fname,figtitle) = read_Valpo(inpath,infile)
    else:
        print('Unrecognized fmt =',fmt)
        df = pd.DataFrame()
        out_fname = ''
        figtitle = ''
    return df, out_fname, figtitle

def read_Albany(inpath,infile,fmt):

    # parse file name
    (base,ext) = os.path.splitext(infile)
    (category,platform,fileDate) = base.split('.')
    #(stn_id_plot,location) = product.split('-')

    # get file time
    file_time = datetime.datetime.strptime(fileDate,'%Y%m%d%H%M')
    #figure_time = 
    
    if fmt == 'Albany':
        stn_id = 'Albany'
        stn_id_plot = 'UAlbany'
        #out_fname = 'upperair.UALB_sonde.{dt}.Albany_NY_skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
        out_fname = 'upperair.UALB.{dt}.Albany_NY_skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    elif fmt == 'Albany_mobile':
        stn_id = 'Albany_Mobile'
        stn_id_plot = 'UAlbany_Mobile'
        #out_fname = 'upperair.UALB_sonde.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
        out_fname = 'upperair.UALB.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))

    stn_id_file = 'UALB'

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

    # remove duplicate columns (keep first)
    df2 = df.loc[:,~df.T.duplicated(keep='first')]

    # get lat and lon
    latStr = df2.iloc[0]['Lat/N']
    deg = int(latStr[0:2])
    min = int(latStr[3:5])
    sec = float(latStr[6:10])
    dir = latStr[-1]
    lat = dirMinSec2deg(deg,min,sec,dir)

    lonStr = df2.iloc[0]['Long/E']
    deg = int(lonStr[0:3])
    min = int(lonStr[4:6])
    sec = float(lonStr[7:11])
    dir = lonStr[-1]
    lon = dirMinSec2deg(deg,min,sec,dir)

    # define figtitle
    if fmt == 'Albany':
        figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id_plot, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)  
    elif fmt == 'Albany_mobile':
        figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id_plot, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)  
    
    for key in df2.keys():
    	# if missing values in data, have to replace with NaNs and then convet to numeric as it
        # originaly is read strings instead of float64; coerce flag nicely takes anything that
        # cannot be converted and puts NaN
        if key != 'UTC_Date' and key != 'UTC_Time' and key != 'AM-PM' and key != 'Lat/N' and key != 'Long/E':
            if df2[key].dtype.name != 'float64':
                df2[key] = pd.to_numeric(df2[key], errors='coerce')

    # convert wind speed from m/s to knots
    df2['speed'] = df2['speed'] * MPS2KTS

    return df2, out_fname, figtitle

def read_TSPOTINT(inpath,infile,fmt,lat,lon):

    # Since lat and lon are in a separate file, hardwire lat and lon
    #if fmt == 'CMICH':
    #    lat = 43.5762
    #    lon = -84.7735
    #elif fmt == 'UNCA':
    #    lat = 35.6196
    #    lon = -82.5674
    lat = float(lat)
    lon = float(lon)

    # Get datetime info from filename
    (category,snd_id,dateTimeStr,ext) = infile.split('.')
    #datetimeObj = datetime.datetime.strptime(dateTimeStr,'%Y%m%d%H%M')

    # Assign filename and figure title
    file_time = datetime.datetime.strptime(dateTimeStr, file_in_dt_fmt_UWYO)
    out_fname = 'upperair.{snd}.{dt}.skewT.png'.format(snd=snd_id, dt=file_time.strftime(file_out_dt_fmt))
    figtitle = '{snd} {dt} sounding ({lati:.3f}, {long:.3f})'.format(snd=snd_id, dt=file_time.strftime(title_dt_fmt), lati=lat, long=lon)

    # Read column headers
    with io.open(inpath+'/'+infile, encoding="cp1252") as myfile:
        head = [next(myfile) for x in range(3)] 
    columnStr = head[0].strip()
    columns = columnStr.split()
    if 'GPM_MSL' in columns:
        columns = [sub.replace('GPM_MSL', 'height') for sub in columns]
    elif 'GPM_AGL' in columns:
        columns = [sub.replace('GPM_AGL', 'height') for sub in columns]
    columns = [sub.replace('Press', 'pressure') for sub in columns]
    columns = [sub.replace('Temp', 'temperature') for sub in columns]
    columns = [sub.replace('RelHum', 'rh') for sub in columns]
    columns = [sub.replace('WDirn', 'direction') for sub in columns]
    columns = [sub.replace('WSpeed', 'speed') for sub in columns]
    
    # Create dataframe
    # read in file after header as DataFrame
    df = pd.read_csv(inpath+'/'+infile, skiprows=3, delimiter='\s+', header=None)
    df.columns = columns
            
    # Calculate dewpoint and add to dataframe
    dew = []
    for i,val in enumerate(df['rh']):
        dew.append(calculate_dewpoint(df['temperature'][i], df['rh'][i]))
    dewpt = np.asarray(dew)
    df.insert(3, "dewpoint", dewpt, True)
        
    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'),
                   how='all').reset_index(drop=True)

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
    #out_fname = 'upperair.CSU_sonde.{dt}.skewT.png'.format(dt=file_time.strftime(file_out_dt_fmt))
    out_fname = 'upperair.CSU.{dt}.skewT.png'.format(dt=file_time.strftime(file_out_dt_fmt))
    figtitle = 'CSU {dt} sounding ({lat:.3f} {lon:.3f})'.format(dt=file_time.strftime(title_dt_fmt), lat=lat_grnd, lon=lon_grnd)

    return df, out_fname, figtitle

def read_MUtxt_ws(inpath,infile):

    # read file header - assumes this format (whitespace between fields, not tabs):
    #%TITLE%
    #XXX   230119/1740
    #Site: lat=41.66758 lon=-73.89545
    # Saved by user: User on 20230119/1826 UTC
    #   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD
    #-------------------------------------------------------------------
    #%RAW%

    with io.open(inpath+'/'+infile, encoding="cp1252") as myfile:
        head = [next(myfile) for x in range(7)] 
    datetimeStr = head[1].strip()
    datetimeStr = re.sub(' +', ' ', datetimeStr)
    (junk,datetimeStr) = datetimeStr.split(' ')
    datetimeObj = datetime.datetime.strptime(datetimeStr,'%y%m%d/%H%M')
    for line in head:
        if 'lat' in line:
            latlonStr = line
            (junk,latStr,lonStr) = latlonStr.split(' ')
            if ',' in latStr:
                latStr = latStr.replace(',','')
            lat = float(latStr.replace('lat=',''))
            lon = float(lonStr.replace('lon=',''))

    # determine out_fname and figtitle
    # Assumes infile naming convention is: upperair.UMILL_windsonde1.202202250600.txt
    parts = infile.split('.')
    snd_id = parts[1]
    file_time_string = parts[2]
        
    file_time = datetime.datetime.strptime(file_time_string, file_in_dt_fmt_UWYO)
    out_fname = 'upperair.{snd}.{dt}.skewT.png'.format(snd=snd_id, dt=file_time.strftime(file_out_dt_fmt))
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = snd_id, dt = datetimeObj.strftime(title_dt_fmt),
                                                                     lati=lat, long=lon)

    # Define missing_val
    missing_val = -9999.00

    # read in file after header as DataFrame
    df = pd.read_csv(inpath+'/'+infile, skiprows=7, header=None)

    # Create column headings
    df.columns = ['pressure','height','temp','dewpt','wdir','wspd']

    # Remove last row ('%END%')
    while df.iloc[-1]['pressure'] == '%END%':
        numRows = len(df.index)
        df.drop([numRows-1],inplace=True)

    # Remove rows containing -9999.00
    df = df[df.temp != missing_val]
    df = df[df.dewpt != missing_val]

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

def read_sharppy(inpath,infile,fmt,lat,lon):

    # Assumes infile naming convention is one of these:
    # upperair.NCSU_sonde.202212271200.txt
    # upperair.Purdue_sonde.202212271200.txt
    
    # Header looks like this:
    # %TITLE%
    #  SITEID   YYMMDD/HHMM LAT,LON             (This line may not contain LAT,LON)
    #  Saved by user: User on 20200220/2320 UTC (This line may be blank)
    #    LEVEL(mb)   HGHT(m MSL)  TEMP(C)   DWPT(C)   WDIR(deg)   WSPD(kts)
    # ---------------------------------------------------------------------

    # read file header
    with open(inpath+'/'+infile) as myfile:
        head = [next(myfile) for x in range(6)]

    info = head[1].strip()
    parts = info.split()
    if len(parts) >= 2:
        stn_id = parts[0]
        datetimeStr = parts[1]
        if len(parts) == 3:
            latLonStr = parts[2]
            (latStr,lonStr) = latLonStr.split(',')
            # turn lat/lon into floats
            lat = float(latStr[:-1])
            lon = float(lonStr[:-1])
        else:
            #if fmt == 'NCSU':
            #    lat = 35.7818
            #    lon = -78.6766
            #elif fmt == 'Purdue':
            #    lat = 40.4302
            #    lon = -86.9154
            #elif fmt == 'UNDws':
            #    lat = 47.3896
            #    lon = -94.6974
            #elif fmt == 'HsinChu':
            #    lat = 24.8279
            #    lon = 121.0142
            #else:  # Unknown location - use 0/0
            #    lat = 0.0
            #    lon = 0.0
            lat = float(lat)
            lon = float(lon)

    else:
        print('No date/time given ... assign Jan 1, 2000')
        datetimeStr = '000101/0000'

    # get file time
    file_time = datetime.datetime.strptime(datetimeStr,'%y%m%d/%H%M')

    # get stn_id - this overrides header info
    if fmt == 'NCSU':
        stn_id = 'NCSU'
    elif fmt == 'Purdue':
        stn_id = 'PurdueU'
    elif fmt == 'UNDws':
        stn_id = 'UND'
    elif fmt == 'HsinChu':
        stn_id = 'HsinChu_46757'
                
    # define out_fname and figtitle
    if 'ws' in fmt:
        out_fname = 'upperair.{stn}_windsonde.{dt}.skewT.png'.format(stn = stn_id,
                                                                     dt = file_time.strftime(file_out_dt_fmt))
        figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id+'_windsonde',
                                                                         dt = file_time.strftime(title_dt_fmt),
                                                                         lati=lat, long=lon)  
    else:
        #out_fname = 'upperair.{stn}_sonde.{dt}.skewT.png'.format(stn = stn_id,
        #                                                         dt = file_time.strftime(file_out_dt_fmt))
        out_fname = 'upperair.{stn}.{dt}.skewT.png'.format(stn = stn_id,
                                                                 dt = file_time.strftime(file_out_dt_fmt))
        #figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id+'_sonde',
        #                                                                 dt = file_time.strftime(title_dt_fmt),
        #                                                                 lati=lat, long=lon)  
        figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id,
                                                                         dt = file_time.strftime(title_dt_fmt),
                                                                         lati=lat, long=lon)  
    # Create dataframe and rename columns
    df = pd.read_csv(inpath+'/'+infile, skiprows=6)
    df.columns = ['pressure','height','temperature','dewpoint','direction','speed']

    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'),
                   how='all').reset_index(drop=True)

    # Handle missing values in data
    for key in df.keys():
	# Replace missing values with NaNs and then convert to numeric as they are originally read 
        # as strings instead of float64; coerce flag nicely takes anything that cannot be converted
        # and substitutes NaN
        if df[key].dtype.name != 'float64':
            df[key] = pd.to_numeric(df[key], errors='coerce')

    return df, out_fname, figtitle

def read_sharppy2(inpath,infile):

    #### NEED TO ADD LAT AND LON TO FIGURE TITLE ####
    
    # determine out_fname and figtitle
    # Assumes infile naming convention is: OSW20221227_12Z_sharppy.txt
    # Assumes infile naming convention is: upperair.OSW_sonde.202212271200.txt
    (base,ext) = os.path.splitext(infile)
    (category,platform,datetimeStr) = base.split('.')
    snd_id = platform
    (snd_id2,junk) = platform.split('_')
    file_time_string = datetimeStr
        
    file_time = datetime.datetime.strptime(file_time_string, file_in_dt_fmt_UWYO)
    out_fname = 'upperair.{snd}.{dt}.skewT.png'.format(snd=snd_id, dt=file_time.strftime(file_out_dt_fmt))
    figtitle = '{snd} {dt} sounding'.format(snd=snd_id2, dt=file_time.strftime(title_dt_fmt))

    # read in file after header as DataFrame
    df = pd.read_csv(inpath+'/'+infile, skiprows=6, header=None)

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
        #out_fname = 'upperair.SBU_sonde.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
        out_fname = 'upperair.SBU.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    elif fmt == 'SBUnc':
        stn_id = 'SBU'
        #out_fname = 'upperair.SBU_sonde.{dt}.Stonybrook_NY_skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
        out_fname = 'upperair.SBU.{dt}.Stonybrook_NY_skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    elif fmt == 'MUnc':
        stn_id = 'UMILL'
        #out_fname = 'upperair.UMILL_sonde.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
        out_fname = 'upperair.UMILL.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    
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

    # assume filename uses this format:
    # upperair.UILL_sonde.YYYYMMDDhhmm.nc
    (junk1,junk2,fileDate,ext) = infile.split('.')
    
    # determine out_fname and figtitle
    attrs = ds.attrs
    stn_id = 'UILL_Mobile'
    figure_time = datetime.datetime.strptime(attrs['start_datetime'],'%Y-%m-%dT%H:%M:%SZ')
    file_time = datetime.datetime.strptime(fileDate,'%Y%m%d%H%M')
    location_str = attrs['location']
    location_list = location_str.split()
    lat = float(location_list[0])
    lon = float(location_list[3])
    if location_list[-1] == 'west':
        lon *= -1
    if location_list[2] == 'south':
        lat *= -1
    #out_fname = 'upperair.UILL_sonde.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    out_fname = 'upperair.UILL.{dt}.skewT.png'.format(dt = file_time.strftime(file_out_dt_fmt))
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id, dt = figure_time.strftime(title_dt_fmt), lati=lat, long=lon)

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

def read_gfk(inpath,infile):

    # Expecting this header:
    
    # Extended NOAA/GMD preliminary data 06-12-2022 21:17:42 [GMT], SkySonde Client version 1.3.1.0
    # Software written by Allen Jordan and Jim Wendell, NOAA
    #                Header lines = 19
    #                Data columns = 15
    #               Flight number = gfk411
    #                  Date [GMT] = 06-12-2022
    #                  Time [GMT] = 21:35:06
    #             Instrument type = Vaisala Sonde
    #              Vaisala number = 52670
    #                  A/D System = V7
    #
    #
    # THE DATA CONTAINED IN THIS FILE ARE PRELIMINARY
    #  AND SUBJECT TO REPROCESSING AND VERIFICATION
    #  FOR SOUNDING DETAILS REFER TO FILE gfk411.de1
    #
    #
    #  Time,     Press,       Alt,      Temp,      Traw,     Theta,        RH,     TFp V,    TVaisI,  Fl,   GPS lat,   GPS lon,   GPS alt,      Wind,  Wind Dir
    # [min],     [hpa],      [km],   [deg C],   [deg C],       [K],       [%],   [deg C],    [degC],  [],     [deg],     [deg],      [km],     [m/s],     [deg]

    missing_value = 99999.
    
    # read file header
    with io.open(inpath+'/'+infile, encoding="cp1252") as myfile:
        head = [next(myfile) for x in range(19)]
        
    dateStr = head[5].strip()
    timeStr = head[6].strip()
    columns = head[17].strip().split(',')
    for idx in range(0,len(columns)):
        columns[idx] = columns[idx].lstrip()        

    # get file time
    parts = dateStr.split()
    date = parts[-1]
    parts = timeStr.split()
    time = parts[-1]
    datetimeStr = date+' '+time
    file_time = datetime.datetime.strptime(datetimeStr,'%d-%m-%Y %H:%M:%S')

    # get station id
    stn_id_file = 'UND'
    stn_id_plot = 'UND'
        
    # read in data as dataframe object; rename/remove columns as needed
    df = pd.read_csv(inpath+'/'+infile, skiprows=19, encoding="latin1", header=None)
    df.columns = columns
    lat = df['GPS lat'].loc[0]
    lon = df['GPS lon'].loc[0]
    df = df.drop(['Time','Traw','Theta','TFp V','TVaisI','Fl','GPS lat','GPS lon','GPS alt'],axis=1)
    df = df.rename(columns={'Press':'pressure','Alt':'height','Temp':'temperature','RH':'rh','Wind':'speed','Wind Dir':'direction'})

    # remove data after balloon stops to rise
    max_ht = df['height'].max()
    idx_list = df.index[df['height']==max_ht].tolist()
    num_rows_to_drop = len(df) - idx_list[0]
    df.drop(df.tail(num_rows_to_drop).index,inplace = True)
    
    # add dewpoint to dataframe
    dew = []
    for i,val in enumerate(df['rh']):
        dew.append(calculate_dewpoint(df['temperature'][i], df['rh'][i]))
    dewpt = np.asarray(dew)
    df.insert(3, "dewpoint", dewpt, True)

    # replace missing values with NaN's
    df = df.replace(missing_value, np.nan)
    
    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'),
                   how='all').reset_index(drop=True)

    # convert wind speed from m/s to knots
    df['speed'] = df['speed'] * MPS2KTS

    # convert height from km to meters
    df['height'] = df['height'] * KM2M

    # define out_fname and figtitle
    #out_fname = 'upperair.{stn}_sonde.{dt}.skewT.png'.format(stn = stn_id_file, dt = file_time.strftime(file_out_dt_fmt))
    out_fname = 'upperair.{stn}.{dt}.skewT.png'.format(stn = stn_id_file, dt = file_time.strftime(file_out_dt_fmt))
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id_plot, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)

    return df, out_fname, figtitle

def read_gfk_csv(inpath,infile):

    # Expecting this header:
    
    #NOAA/GMD balloon flight data file
    #
    #SkySonde Client Version:, 1.3.1.0
    #Balloon Flight Name:, gfk411
    #Attached Instruments:, Radiosonde
    #Station Name:, Grand Forks, ND
    #Station Latitude:, 47.92131
    #Station Longitude:, -97.08764
    #Station Altitude:, 256.00 [m]
    #Station Pressure:, 99999 [mb]
    #Station Temperature:, 99999 [deg C]
    #Station Humidity:, 99999 [% RH]
    #Radiosonde Serial Number / ID:, 69116
    #Pressure/Altitude Source:, iMet Radiosonde
    #Solar Radiation Correction:, Yes
    #
    #date [y-m-d GMT], time [h:m:s GMT], milliseconds, seconds since midnight [GMT], elapsed minutes, altitude (from iMet PTU) [km], iMet pressure [mb], iMet air temperature (corrected) [deg C], iMet air temperature (raw) [deg C], iMet humidity [RH %], iMet frostpoint [deg C], iMet internal temperature [deg C], iMet battery voltage [V], iMet theta [K], iMet temperature (of pressure sensor) [deg C], iMet temperature (of humidity sensor) [deg C], iMet ascent rate [m/s], iMet water vapor mixing ratio [ppmv], iMet total column water [mm], GPS latitude, GPS longitude, GPS altitude [km], GPS num satellites, GPS pressure [mb], GPS wind speed [m/s], GPS wind direction [deg], GPS ascent rate [m/s], GPS(X) east velocity [m/s], GPS(X) north velocity [m/s], GPS(X) up velocity [m/s], GPS time [h:m:s GMT], GPS heading from launch [deg], GPS elevation angle from launch [deg], GPS distance from launch [km], predicted landing latitude, predicted landing longitude, predicted time to landing [min]

    missing_value = 99999.
    
    # read file header
    with io.open(inpath+'/'+infile, encoding="cp1252") as myfile:
        head = [next(myfile) for x in range(16)]
        
    latStr = head[6].strip()
    lonStr = head[7].strip()
    lat = float(latStr.split()[-1])
    lon = float(lonStr.split()[-1])

    # get station id
    stn_id_file = 'UND'
    stn_id_plot = 'UND'
        
    # read in data as dataframe object; rename/remove columns as needed
    df = pd.read_csv(inpath+'/'+infile, skiprows=16, encoding="latin1")

    # get file time
    date = df['date [y-m-d GMT]'][0]
    time = df[' time [h:m:s GMT]'][0]
    datetimeStr = date+' '+time
    file_time = datetime.datetime.strptime(datetimeStr,'%Y-%m-%d %H:%M:%S')

    df = df.drop(['date [y-m-d GMT]',' time [h:m:s GMT]',' milliseconds',' seconds since midnight [GMT]',' elapsed minutes',
                  ' iMet air temperature (raw) [deg C]',' iMet frostpoint [deg C]',' iMet internal temperature [deg C]',
                  ' iMet battery voltage [V]',' iMet theta [K]',' iMet temperature (of pressure sensor) [deg C]',
                  ' iMet temperature (of humidity sensor) [deg C]',' iMet ascent rate [m/s]',' iMet water vapor mixing ratio [ppmv]',
                  ' iMet total column water [mm]', ' GPS latitude', ' GPS longitude',' GPS altitude [km]',' GPS num satellites',
                  ' GPS pressure [mb]',' GPS ascent rate [m/s]', ' GPS(X) east velocity [m/s]',' GPS(X) north velocity [m/s]',
                  ' GPS(X) up velocity [m/s]',' GPS time [h:m:s GMT]',' GPS heading from launch [deg]',
                  ' GPS elevation angle from launch [deg]',' GPS distance from launch [km]',' predicted landing latitude',
                  ' predicted landing longitude',' predicted time to landing [min]'],axis=1)
    df = df.rename(columns={' iMet pressure [mb]':'pressure',' altitude (from iMet PTU) [km]':'height',
                            ' iMet air temperature (corrected) [deg C]':'temperature',' iMet humidity [RH %]':'rh',
                            ' GPS wind speed [m/s]':'speed',' GPS wind direction [deg]':'direction'})

    # add dewpoint to dataframe
    dew = []
    for i,val in enumerate(df['rh']):
        dew.append(calculate_dewpoint(df['temperature'][i], df['rh'][i]))
    dewpt = np.asarray(dew)
    df.insert(3, "dewpoint", dewpt, True)

    # replace missing values with NaN's
    df = df.replace(missing_value, np.nan)
    
    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'),
                   how='all').reset_index(drop=True)

    # convert wind speed from m/s to knots
    df['speed'] = df['speed'] * MPS2KTS

    # convert height from km to meters
    df['height'] = df['height'] * KM2M

    # define out_fname and figtitle
    #out_fname = 'upperair.{stn}_sonde.{dt}.skewT.png'.format(stn = stn_id_file, dt = file_time.strftime(file_out_dt_fmt))
    out_fname = 'upperair.{stn}.{dt}.skewT.png'.format(stn = stn_id_file, dt = file_time.strftime(file_out_dt_fmt))
    figtitle = '{stn} {dt} sounding ({lati:.3f}, {long:.3f})'.format(stn = stn_id_plot, dt = file_time.strftime(title_dt_fmt), lati=lat, long=lon)

    return df, out_fname, figtitle

def read_UWYO(inpath,infile):

    # Get out_fname and figtitle
    parts = infile.split('.')
    file_time_string = parts[2]
    snd_id = parts[3]
    snd_id = siteList[snd_id]
    
    file_time = datetime.datetime.strptime(file_time_string, file_in_dt_fmt_UWYO)
    out_fname = 'upperair.SkewT.{dt}.{snd}.png'.format(dt=file_time.strftime(file_out_dt_fmt),snd=snd_id)
    figtitle = '{snd} {dt} sounding'.format(snd=snd_id, dt=file_time.strftime(title_dt_fmt))
    
    # Read infile into array of lines
    with open(inpath+'/'+infile) as fid:
        lines = [line for line in fid.readlines() if line.strip()]
        nlines = len(lines)

    # Get required info from header
    header = lines[1]
    fields = lines[3].split()
    nfields = len(fields)
    units = lines[4].split()

    # Create temp file with only data
    tmpFile = '/tmp/'+infile+'.tmp'
    numRecs = 0
    fnew = open(tmpFile,'w')
    lines_without_header = lines[6:]
    for idx,iline in enumerate(lines_without_header,0):
        if not iline.startswith('Station'):
            parts = re.findall(r'\S+', iline)
            if len(parts) == nfields:
                fnew.write(iline)
                numRecs = numRecs + 1
        else:
            break
    fnew.close()

    if numRecs > 1:
    
        # Read data from tmpFile
        if header[:5] == '00000':
            # WRF profile
            station = '-99999'
            Longitude = float(header.split()[5].strip(","))
            Latitude = float(header.split()[6])
            sounding_date = header.split()[-1]
        else:
            # UWYO profile
            station = header[:5]
            dstr = (' ').join(header.split()[-4:])
            sounding_date = datetime.datetime.strptime(dstr, "%HZ %d %b %Y").strftime("%Y-%m-%d_%H:%M:%S") 

        #file_data = np.genfromtxt(inpath+'/'+infile, skip_header=16, skip_footer=nlines_footer)        
        file_data = np.genfromtxt(tmpFile)
        os.remove(tmpFile)
        
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

    else:
        os.remove(tmpFile)
        df = pd.DataFrame()
        
    return df, out_fname, figtitle

def read_Valpo(inpath,infile):

    # Get out_fname and figtitle
    latFound = False
    lonFound = False
    with open(inpath+'/'+infile, 'r') as fp:
        while not latFound or not lonFound:
            line = fp.readline().strip()
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



