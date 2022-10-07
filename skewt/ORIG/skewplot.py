# Brody Fuchs, CSU, November 2015
# brfuchs@atmos.colostate.edu


## Update the skewT code to be smarter using classes and what not
# use the KUILsounding.txt as an example or KBMXsounding.txt

# Start out by just plotting the raw data first

import numpy as np
import matplotlib.pyplot as plt
from skewPy import SkewT
import os
from datetime import datetime
import glob
import sys
import re
import ntpath

##### PARAMETERS YOU CAN CHANGE HERE #########

file_path = '/home/disk/funnel/olympex/archive/research/text_sounding' # location of sounding files
out_path = '/home/disk/funnel/olympex/archive/research/skewt' # where you want to output skewTs
prefix = 'research' # the first part of the sounding raw data file, need it to search for possible files
                    # other possibilities are 'K' for K*** UWYO soundings, could also be 'EDT' for Vaisala soundings
fmt = 'EC' # This can currently be 'EC' for Canadian format, 'UWYO' for Wyoming soundings or 'EDT' for Vaisala CSU soundings
station_name = None # put this in if want to override the station name in the title of the skewT, default is fmt above
                                    # otherwise just put None

########### DON'T NEED TO WORRY ABOUT ANYTHING BELOW HERE ############


if len(sys.argv) > 1:
    if 'debug' in sys.argv:
        debug = True
    else: 
        debug = False

    if 'number' in sys.argv:
        number_loc = sys.argv.index('number')
        print sys.argv[number_loc + 1]
        sounding_files = [glob.glob('%s/%s*'%(file_path, prefix))[int(sys.argv[number_loc+1])]]
        print '\nonly processing %s in %s'%(sounding_files[0], file_path)

    elif 'last' in sys.argv: 
        sounding_files = [glob.glob('%s/%s*'%(file_path, prefix))[-1]]
        print '\nonly processing most recent file in %s\n\nIf you want to process all files in directory, just type: python skew.py\n'%file_path
    else: 
        sounding_files = glob.glob('%s/%s*'%(file_path, prefix))
        print '\nProcessing all files in %s \n\nIf you only want to process the most recent file, just type: python skew.py last\n'%file_path

else: 
    sounding_files = sorted(glob.glob('%s/%s*'%(file_path, prefix)))
    print sounding_files
    debug = False
    print '\nProcessing all files in %s \n\nIf you only want to process the most recent file, just type: python skew.py last\n'%file_path


for fname in sounding_files:
    #file_title = os.path.basename(fname)[4:-4]

    print 'Processing %s'%os.path.basename(fname)

    base_fname = ntpath.basename(fname)
    
    fnameParts = re.split('\.',base_fname)
    datetime = fnameParts[2]
    date = datetime[:8]
    print 'date = %s\n'%date
    try:
        os.stat('%s/%s'%(out_path,date))
    except:
        print 'Create new directory for %s\n'%date
        os.mkdir('%s/%s'%(out_path,date))

    S = SkewT.Sounding(fname, fmt = fmt, station_name = station_name)
    S.plot_skewt(parcel = True, parcel_draw = True)

    #plt.savefig('%s/%s_sounding.png'%(out_path, S.sounding_date), dpi = 150)
    plt.savefig('%s/%s/research.skewt.%s.ec.png'%(out_path, date, S.sounding_date), dpi = 150)
    plt.close('all')

print '\n%d total soundings in %s\n'%(len(glob.glob('%s/%s*'%((file_path, prefix)))), file_path)


if debug: # this code is to do some testing before going into the SkewT module


    ptop = 150
    pbot = 400

    trop_weights = {'theta': 1.0, 'dtdp': 1.0, 'T': 1.0, 'dthetadp': 1.0, 'drhdp': 0.6, 'pres': 1.0}

    smooth_temp = SkewT.running_mean(S.data['temp'], 3, ntimes = 30)

    dp = np.gradient(S.data['pres'])

    #lapse = 1000.0*np.diff(smooth)/np.diff(S.data['hght'])

    deriv = np.gradient(smooth_temp, dp)
    deriv_smooth = SkewT.running_mean(deriv, 3, ntimes = 30)


    deriv2 = np.gradient(deriv_smooth, dp)
    deriv2_smooth = 100.0*SkewT.running_mean(deriv2, 5, ntimes = 20)

    ti = SkewT.tropopause_index(S.data['temp'])
    pt = S.data['pres'][ti]


    fig, ax = plt.subplots(1, 3, figsize = (15, 6))

    ax[0].plot(S.data['temp'], S.data['pres'], 'red', label = 'raw data')
    #plt.plot(smooth, S.data['pres'][bad_first:]-delay, 'black', label = 'smoothed')
    #plt.plot(deriv*100, S.data['pres'][bad_first:]-delay, 'blue', label = 'derivative')
    ax[0].plot(smooth_temp, S.data['pres'], 'black', label = 'smoothed')
    #ax[0].plot(100*deriv, S.data['pres'], 'blue', label = 'derivative')
    axderiv = ax[0].twiny()
    axderiv.plot(deriv_smooth, S.data['pres'], 'blue', label = 'smooth derivative')
    #ax[0].plot(deriv2, S.data['pres'], 'green', label = '2nd derivative')
    #ax[0].plot(deriv2_smooth, S.data['pres'], 'green', label = 'smooth 2nd der')
    #plt.plot(lapse, S.data['pres'][:-1], 'brown', label = 'lapse rate')

    #ax[0].axhline(y = pt, color = 'black')
    #ax[0].axvline(x = 0, color = 'black')
    ax[0].set_ylim(pbot, ptop)
    ax[0].set_xlim(-80, 30)
    ax[0].legend(loc = 'best')
    ax[0].grid(True)
    axderiv.set_xlim(-1,1)
    axderiv.axvline(x = 0, color = 'blue', linestyle = 'dashed')

    theta = SkewT.running_mean(SkewT.Theta(273.15+S.data['temp'], S.data['pres']*100.0), 3, ntimes = 50)

    dtheta = -1.0*SkewT.running_mean(np.gradient(theta, dp), 3, ntimes = 10)
    ddtheta = SkewT.running_mean(np.gradient(dtheta, dp), 3, ntimes = 30)
    axd = ax[1].twiny()

    ax[1].plot(theta, S.data['pres'], color = 'black')
    ax[1].set_ylim(pbot, ptop)
    ax[1].set_xlim(270, 450)
    #ax[1].axhline(y = 320)
    axd.plot(dtheta, S.data['pres'], 'blue', label = 'dtheta')
    #axd.plot(ddtheta, S.data['pres'], 'brown', label = '2nd deriv theta')
    #axd.plot(dtheta**3/1000, S.data['pres'], 'green', label = 'deriv cubed')
    
    axd.set_xlim(-1, 3)
    axd.axvline(x = 0, color = 'blue', linestyle = 'dashed')
    axd.legend(loc = 'lower right')


    # last one will be a humidity one
    rh = SkewT.running_mean(SkewT.RH(S.data['temp'], S.data['dwpt']), 3, ntimes = 30)
    drh = -1.0*np.gradient(rh, dp)

    ax[2].plot(rh, S.data['pres'])

    axdh = ax[2].twiny()

    axdh.plot(drh, S.data['pres'], 'green')


    ax[2].set_ylim(pbot, ptop)
    axdh.set_xlim(-2, 2)

    input_data = {'theta': theta, 'dtdp': deriv_smooth, 'T': smooth_temp, 'dthetadp': dtheta, 'drhdp': drh, 'pres': S.data['pres']}
    trop_pressure = SkewT.fuzzy_tropopause(input_data, trop_weights)

    trop_score = SkewT.score_calc(input_data, weights = trop_weights)

    axd.plot(trop_score, S.data['pres'], color = 'darkorange', label = 'tropopause score')

    pressure_fix = 0

    ax[0].axhline(y = trop_pressure)
    ax[0].axhline(y = trop_pressure - pressure_fix, color = 'green')
    ax[1].axhline(y = trop_pressure)
    ax[1].axhline(y = trop_pressure - pressure_fix, color = 'green')
    ax[2].axhline(y = trop_pressure)
    ax[2].axhline(y = trop_pressure - pressure_fix, color = 'green')

    print 'tropopause score maximized at: %d mb, adjusted: %d mb'%(trop_pressure, trop_pressure - pressure_fix)


    plt.tight_layout()
    plt.show()


