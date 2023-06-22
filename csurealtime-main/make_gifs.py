import numpy as np 
import os
import glob

# this script just makes animated gifs of various things

base_path = '/home/rmet/SPURS2/realtime_radar'

nimages = 12
rain_fig_check = sorted(glob.glob('%s/figures/rain/*.png'%(base_path)))
ppi_fig_check = sorted(glob.glob('%s/figures/raw_ppi/*.png'%(base_path)))
large_ppi_fig_check = sorted(glob.glob('%s/figures/large_dbz_vel_ppi/*.png'%(base_path)))

print 'making rain gif'
if len(rain_fig_check): # means there are files there
	if len(rain_fig_check) > nimages:
	    valid_rain_files = rain_fig_check[-1*nimages:]
	else:
	    valid_rain_files = deepcopy(rain_fig_check)

	# now need to make a string out of it
	valid_rain_file_string = ' '.join(valid_rain_files)
	#print 'valid rain file string: {}'.format(valid_rain_file_string)
	os.system('convert -delay 60 -loop 0 %s %s/figures/animated_gifs/rain_latest.gif'%(valid_rain_file_string, base_path))


print 'making ppi gif'
if len(ppi_fig_check): # means there are files there
    if len(ppi_fig_check) > nimages:
        valid_ppi_files = ppi_fig_check[-1*nimages:]
    else:
        valid_ppi_files = deepcopy(ppi_fig_check)

    # now need to make a string out of it
    valid_ppi_file_string = ' '.join(valid_ppi_files)
    #print 'valid rain file string: {}'.format(valid_rain_file_string)
    os.system('convert -delay 60 -loop 0 %s %s/figures/animated_gifs/ppi_latest.gif'%(valid_ppi_file_string, base_path))

print 'making large ppi gif'
if len(large_ppi_fig_check): # means there are files there
    if len(large_ppi_fig_check) > nimages:
        valid_large_ppi_files = large_ppi_fig_check[-1*nimages:]
    else:
        valid_large_ppi_files = deepcopy(large_ppi_fig_check)

    # now need to make a string out of it
    valid_large_ppi_file_string = ' '.join(valid_large_ppi_files)
    #print 'valid rain file string: {}'.format(valid_rain_file_string)
    os.system('convert -delay 60 -loop 0 %s %s/figures/animated_gifs/large_ppi_latest.gif'%(valid_large_ppi_file_string, base_path))