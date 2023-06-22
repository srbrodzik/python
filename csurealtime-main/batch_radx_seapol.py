# Brody Fuchs
# just hacking the batch shell scripting here

import numpy as np
import glob
import os
import argparse
import glob
import sys
import subprocess as sp
import time
import Config
import gentools as gtools 
import plot_tools as ptools 
import radar_tools as radtools 

# This is code that is effectively a shell script that calls RadXConvert which grids radar data
# this code is called after we've QCed the data with qc_radar.py

try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception:
    base_path = os.path.dirname(os.path.realpath('__file__'))


# parse the command line arguments
parser = argparse.ArgumentParser(description='Put in a file to be processed')


parser.add_argument('--file', action="store", dest="file")
parser.add_argument('--lastfile', action="store", dest="lastfile", type=bool, default=True)
parser.add_argument('--allfiles', action="store", dest="allfiles", type=bool, default=False)
parser.add_argument('--realtime', action="store", dest="realtime", type=bool, default=False)
parser.add_argument('--config', action="store", dest="config", default=None)


pargs = parser.parse_args()

# Parameters


if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v


outpath = ycfg['grid_data_path']
#inpath = '%s/%s/ppi'%(base_path, ycfg['path2'])
inpath = '%s/%s/ppi'%(base_path, 'realtime_drops')
params_file = '%s/Radx2Grid.cpol.relampago.cart'%(base_path)
radx_cmd = '/home/radaruser/lrose/bin/Radx2Grid'

out_prefix = ycfg['grid_prefix']
in_prefix = ycfg['qc_prefix']
file_type = 'nc'

print 'out prefix: {}'.format(out_prefix)

mode = 'last'

if pargs.realtime:
	print 'RadX in realtime'
	pass
elif pargs.file is not None:
	mode = 'file'
	files = [pargs.file]

elif pargs.allfiles:
	mode = 'allfiles'
	files = glob.glob('%s/*.%s'%(inpath, file_type))

elif pargs.lastfile:
	mode = 'file'
	files = [glob.glob('%s/*.%s'%(inpath, file_type))[-1]]

else:
	print 'No files specified or found, exiting.'
	sys.exit()




print files

for f in files:

	# the idea here is to check if this file has already been converted, if so then don't convert it again
	fbase = os.path.basename(f)
	stripped_fbase = fbase[len(in_prefix):-1*(len(file_type)+1)]

	print 'input file base: {}'.format(stripped_fbase)

	stripped_date = stripped_fbase[:8]
	print 'stripped date: {}'.format(stripped_date)

	# ***** THIS IS EXTREMELY ANNOYING ********
	# The datetime string gets CHANGED by RadXConvert for some dumbass reason
	# So I grab the datetime string from the file before it gets converted...

#	if len(outpath_file_check) == 0:
	if True:
		# We first call RadX2Grid
		radgrid_cmd = '%s -params %s -f %s -outdir %s -name_start'%(radx_cmd, params_file, f, outpath)
		os.system(radgrid_cmd)

		# Since radx has some weird naming stuff, we're gonna move the file
		# use subprocess to see the most recent file in the output directory
		# This subprocess stuff allows me to grab the name of the most recent file in the outpath directory
		# which should be the file we just converted

		radx_old_fname = os.path.basename(gtools.recent_file_glob('%s/%s'%(outpath, stripped_date)))

		print 'radx out filename: {}'.format(radx_old_fname)

		# Here we grab the name of the RadX output file so we can replace it with the name before it got messed up
		# during conversion
		prefix_ind = radx_old_fname.index(out_prefix) + len(out_prefix)
		radx_fbase = radx_old_fname[prefix_ind:-1*(len(file_type)+1)]

		stripped_radxdate = radx_fbase[:8]
		radx_new_fname = radx_old_fname.replace(radx_fbase, stripped_fbase)


		# then need to move the file using os.system
		os.system('mv %s/%s/%s %s/%s'%(outpath, stripped_radxdate, radx_old_fname, outpath, radx_new_fname))
		latest_data_check = glob.glob('%s/*latest_data_info*'%(outpath))
		if len(latest_data_check): # I'm not sure what these files are but I don't want them
								# so let's blow them away
			print 'Found some latest_data files, deleting them'
			os.system('rm %s/*latest_data_info*'%(outpath))

		time.sleep(1)


	else:
		print 'file: %s has already been gridded'%(stripped_fb)











