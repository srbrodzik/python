# -*- coding: utf-8 -*-
"""
Submit this to loop through all chivo PPI files
"""

#############
from chivo_pf_id_v1_maxdbz import PF_finder_max
import os
import fnmatch
#############
from netCDF4 import Dataset
import datetime
import numpy as np
from hyu_convstrat import conv_stra_sep
from hid_iwp_chivo import get_csu_hid

def read_files(rad_file):
    
    data = Dataset(rad_file, 'r')

    epoch_time = np.array(data.variables['time'])[0]
    datetime_time = datetime.datetime.fromtimestamp(epoch_time)
    hours_added = datetime.timedelta(hours = 7)
    date = (datetime_time + hours_added)
    date = date.strftime("%Y%m%d_%H%M%S")
    print (date)

    x = np.array(data.variables['x0'])
    y = np.array(data.variables['y0'])
    z = np.array(data.variables['z0'])
    reflect = np.array(np.squeeze(data.variables['DZ_qc']))
    #hydro = np.array(np.squeeze(data.variables['HID_qc']))
    hydro, iwp, mw, mi = get_csu_hid(rad_file, snd_times, filenames)
    
    cs_all = []
    x1,y1 = np.meshgrid(x,y)

    #want the max conv/strat between 3 and 5 km (corresponds to z[5:10])
    for z1 in np.arange(5, 10, 1):
        cs,cc,bkgnd = conv_stra_sep(reflect[z1,:,:], y1, x1, 1, 'C', 'SHY')

        cs[np.where(reflect[z1,:,:] <= 0)] = -1
        cs[np.where(np.isnan(reflect[z1,:,:]))] = -1

        cs_all.append(cs)
        
    csmax = np.max(cs_all, axis=0)
    
    return date, x, y, z, reflect, hydro, iwp, csmax
#############
import pandas as pd 

def write_output(out, date):

    #length, width for out[0], out[20], out[27], out[34], out[40]
    length = []
    width = []

    for g in np.arange(0, len(out[0]), 1):
        el = out[0][g]
        len1 = el.width
        wid1 = el.height
        if (len1 > wid1):
            length.append(len1)
            width.append(wid1)
        else:
            width.append(len1)
            length.append(wid1)

    dcc_length = np.zeros((len(out[0])), dtype='object')
    dcc_width = np.zeros((len(out[0])), dtype='object')

    for g in np.arange(0, len(out[0]), 1):
        elc = out[20][g]                  # change 20 to the correct number (DCC_ellipse) if changing output
        if not elc:
            clen = 0
            cwid = 0  
            dcc_length[g] = clen
            dcc_width[g] = cwid
        else:
            dcc_length1 = np.zeros((len(elc)))
            dcc_width1 = np.zeros((len(elc)))
            for m in np.arange(0, len(elc), 1):
                elc_noempty = elc[m]
                clen = elc_noempty.width
                cwid = elc_noempty.height
                if (clen > cwid):
                    dcc_length1[m] = clen
                    dcc_width1[m] = cwid
                else:
                    dcc_width1[m]  = clen
                    dcc_length1[m] = cwid
            dcc_length[g] = dcc_length1
            dcc_width[g] = dcc_width1

    dwcc_length = np.zeros((len(out[0])), dtype='object')
    dwcc_width = np.zeros((len(out[0])), dtype='object')

    for g in np.arange(0, len(out[0]), 1):
        elc = out[27][g]               # change 27 to the correct number (DWCC_ellipse) if changing output
        if not elc:
            clen = 0
            cwid = 0
            dwcc_length[g] = clen
            dwcc_width[g] = cwid
        else:
            dwcc_length1 = np.zeros((len(elc)))
            dwcc_width1 = np.zeros((len(elc)))
            for m in np.arange(0, len(elc), 1):
                elc_noempty = elc[m]
                clen = elc_noempty.width
                cwid = elc_noempty.height
                if (clen > cwid):
                    dwcc_length1[m] = clen
                    dwcc_width1[m] = cwid
                else:
                    dwcc_width1[m]  = clen
                    dwcc_length1[m] = cwid
            dwcc_length[g] = dwcc_length1
            dwcc_width[g] = dwcc_width1

    wcc_length = np.zeros((len(out[0])), dtype='object')
    wcc_width = np.zeros((len(out[0])), dtype='object')

    for g in np.arange(0, len(out[0]), 1):
        elc = out[34][g]                 # change 34 to the correct number (WCC_ellipse) if changing output
        if not elc:
            clen = 0
            cwid = 0
            wcc_length[g] = clen
            wcc_width[g] = cwid
        else:
            wcc_length1 = np.zeros((len(elc)))
            wcc_width1 = np.zeros((len(elc)))
            for m in np.arange(0, len(elc), 1):
                elc_noempty = elc[m]
                clen = elc_noempty.width
                cwid = elc_noempty.height
                if (clen > cwid):
                    wcc_length1[m] = clen
                    wcc_width1[m] = cwid
                else:
                    wcc_width1[m]  = clen
                    wcc_length1[m] = cwid
            wcc_length[g] = wcc_length1
            wcc_width[g] = wcc_width1


    bsr_length = np.zeros((len(out[0])), dtype='object')
    bsr_width = np.zeros((len(out[0])), dtype='object')

    for g in np.arange(0, len(out[0]), 1):
        elc = out[40][g]               # change 40 to the correct number (BSR_ellipse) if changing output
        if not elc:
            clen = 0
            cwid = 0
            bsr_length[g] = clen
            bsr_width[g] = cwid
        else:
            bsr_length1 = np.zeros((len(elc)))
            bsr_width1 = np.zeros((len(elc)))
            for m in np.arange(0, len(elc), 1):
                elc_noempty = elc[m]
                clen = elc_noempty.width
                cwid = elc_noempty.height
                if (clen > cwid):
                    bsr_length1[m] = clen
                    bsr_width1[m] = cwid
                else:
                    bsr_width1[m]  = clen
                    bsr_length1[m] = cwid
            bsr_length[g] = bsr_length1
            bsr_width[g] = bsr_width1

    
    
    #create dictionary of all variables
    dic= ({'xloc':out[1], 'yloc':out[2], 'pf_length':length, 'pf_width':width, 'pf_area':out[3], 'convective_area':out[4],
           'stratiform_area':out[5], 'mean_dbz_by_altitude':out[6], 'max_dbz_by_altitude':out[7], 'echo_top':out[8], 
           '40dbz_top':out[9], 'mean_rainrate':out[10], 'max_rainrate':out[11], 'mean_convective_rainrate':out[12], 
           'mean_stratiform_rainrate':out[13], 'mean_rainrate_mdbz':out[14], 'max_rainrate_mdbz':out[15], 
           'mean_convective_rainrate_mdbz':out[16], 'mean_stratiform_rainrate_mdbz':out[17], 'iwp_max':out[18], 'iwp_mean':out[19],
           'dcc_xloc':out[22], 'dcc_yloc':out[23], 'dcc_length':dcc_length, 'dcc_width':dcc_width,
           'dcc_area':out[21], 'dcc_dbz_at_10km':out[24], 'dcc_mean_rainrate':out[25], 'dcc_mean_rainrate_mdbz':out[26],
           'dwcc_xloc':out[29], 'dwcc_yloc':out[30], 'dwcc_length':dwcc_length, 'dwcc_width':dwcc_width,
           'dwcc_area':out[28], 'dwcc_dbz_at_10km':out[31], 'dwcc_mean_rainrate':out[32],
           'dwcc_mean_rainrate_mdbz':out[33], 'wcc_xloc':out[36], 'wcc_yloc':out[37], 'wcc_length':wcc_length,
           'wcc_width':wcc_width, 'wcc_area':out[35], 'wcc_mean_rainrate':out[38], 'wcc_mean_rainrate_mdbz':out[39],
           'bsr_xloc':out[42], 'bsr_yloc':out[43], 'bsr_length':bsr_length, 'bsr_width':bsr_width, 'bsr_area':out[41],
           'bsr_mean_rainrate':out[44], 'bsr_mean_rainrate_mdbz':out[45], '35dbz_xloc':out[46], '35dbz_yloc':out[47],
           '35dbz_volume':out[48], 'graupel_xloc':out[49], 'graupel_yloc':out[50], 'graupel_volume':out[51],
           'graupel_max_height':out[52], 'percent_graupel':out[53], 'percent_rain':out[54], 'percent_snow':out[55],
           'percent_ice': out[56], 'percent_hail':out[57]})
       
    df= pd.DataFrame(dic)
    
    ind = np.where(df['echo_top'] > 0)
    out1 = df.iloc[ind[0]] 
    out1['time'] = date
    
    #create a csv file from the dictionary
    out1.to_csv("/rasmussen-scratch2/mrocque/research/relampago/chivo_trmm_pfs/final_steiner/csu_chivo_1kmgrid_pfs_ppi_shy_"+date+".csv", 
                columns=['time', 'xloc', 'yloc', 'pf_length', 'pf_width', 'pf_area', 'convective_area', 'stratiform_area',
                         'mean_dbz_by_altitude', 'max_dbz_by_altitude', 'echo_top', '40dbz_top', 'mean_rainrate', 
                         'max_rainrate', 'mean_convective_rainrate', 'mean_stratiform_rainrate', 'mean_rainrate_mdbz',
                         'max_rainrate_mdbz', 'mean_convective_rainrate_mdbz', 'mean_stratiform_rainrate_mdbz', 'iwp_max', 'iwp_mean',
                         'dcc_xloc', 'dcc_yloc', 'dcc_length', 'dcc_width', 'dcc_area', 'dcc_dbz_at_10km', 'dcc_mean_rainrate',
                         'dcc_mean_rainrate_mdbz', 'dwcc_xloc', 'dwcc_yloc', 'dwcc_length', 'dwcc_width', 'dwcc_area',
                         'dwcc_dbz_at_10km', 'dwcc_mean_rainrate', 'dwcc_mean_rainrate_mdbz', 'wcc_xloc', 'wcc_yloc', 'wcc_length',
                         'wcc_width', 'wcc_area', 'wcc_mean_rainrate', 'wcc_mean_rainrate_mdbz', 'bsr_xloc', 'bsr_yloc',
                         'bsr_length', 'bsr_width', 'bsr_area', 'bsr_mean_rainrate', 'bsr_mean_rainrate_mdbz', '35dbz_xloc',
                         '35dbz_yloc', '35dbz_volume', 'graupel_xloc', 'graupel_yloc', 'graupel_volume', 'graupel_max_height',
                         'percent_graupel', 'percent_rain', 'percent_snow', 'percent_ice', 'percent_hail'])

##################  find the current files (useful if you need to restart the code)    
curr_files = []

for path, dirs, files in sorted(os.walk('/rasmussen-scratch2/mrocque/research/relampago/chivo_trmm_pfs/final_steiner/')):
    for file in sorted(files):
        if fnmatch.fnmatch(file, '*.csv'):
            curr_files.append(file)
            
##################  grab all the sounding times
snd_times = []
filenames = []
for filename in os.listdir('/rasmussen-scratch/mrocque/research/relampago/CACTI_ARM_soundings/netcdf/'):
    filenames.append(filename)
    str2 = (filename[18:33])
    snd_times.append(datetime.datetime.strptime(str2, '%Y%m%d.%H%M%S'))
    
#put them in order
filenames.sort()
snd_times.sort()

############# loop thru a directory and compare the files to current files, calculate PFs, write output to csv file
for path,dirs,files in sorted(os.walk('/rasmussen-scratch/krasmussen/DATA/RELAMPAGO/CHIVO_interpolated_data_fromStacy/ppi/')):
        for file in sorted(files):
            if fnmatch.fnmatch(file,'*.nc'):
                fullname = os.path.join(path,file)
                if 'csu_chivo_1kmgrid_pfs_ppi_'+str(file[20:35])+'.csv' not in curr_files:
                    print (file[20:35])
                    date, x, y, z, reflect, hydro, iwp, csmax = read_files(fullname)
                    out = PF_finder_max(reflect, x, y, z, csmax, hydro, iwp, 3.0, 0, 40, 10, 1000, 10000)
                    write_output(out, date)
                    print ('done with ' + date)  
# ##############################################################################