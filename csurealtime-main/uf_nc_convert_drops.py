# Brody Fuchs
# just hacking the batch shell scripting here
# This will take output from the original QC'ed data and convert it to UF, then run
# DROPS on that uf file, then convert that uf file back to cfradial for subsequent use in the processing pipeline

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


drops_config_file = ycfg['drops_config_file']


inpath = '%s/%s/ppi'%(base_path, ycfg['qc_data_path'])

outpath = '%s/%s'%(base_path, ycfg['path3'])

params_file = '%s/Radx2Grid.seapol.spurs2.cart'%(base_path)
out_prefix = ''
in_prefix = ycfg['prefix2']
file_type = 'nc'


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


	print 'f: {}'.format(f)

	print 'input filename: {}'.format(stripped_fbase)

	stripped_date = stripped_fbase[:8]
	stripped_month = stripped_fbase[:6]
	print 'stripped date: {}'.format(stripped_date)

	# ***** THIS IS EXTREMELY ANNOYING ********
	# The datetime string gets CHANGED by RadXConvert for some dumbass reason
	# So I grab the datetime string from the file before it gets converted...

#	if len(outpath_file_check) == 0:
	if True:
		# We first call RadX2Grid
		#os.system('Radx2Grid -params %s -f %s -outdir %s -name_start'%(params_file, f, outpath))
		radx_cmd = 'RadxConvert -f %s -uf -outdir %s/realtime_drops'%(f, base_path)
		
		os.system(radx_cmd)

		print 'DONE WITH FIRST CONVERSION TO UF'

		# that puts the file in the realtime_drops/$stripped_date directory
		# we need to run DROPS on that file, but where do we put it??

		blah_path = '%s/realtime_drops/%s'%(base_path, stripped_date)
		#print 'blah path: {}'.format(blah_path)

		#radx_old_fname = gtools.recent_file(blah_path)
		radx_old_fname = os.path.basename(gtools.recent_file_glob(blah_path, index=-1, before_wildcard='*%s'%(stripped_month)))

		#print 'radx out filename: {}'.format(radx_old_fname)


		# Here we grab the name of the RadX output file so we can replace it with the name before it got messed up
		# during conversion
		prefix_ind = radx_old_fname.index(out_prefix) + len(out_prefix)
		radx_fbase = radx_old_fname[prefix_ind:-1*(len(file_type)+1)]


		stripped_radxdate = radx_fbase[:8]
		radx_new_fname = radx_old_fname.replace(radx_fbase, stripped_fbase)
		#radx_new_fname = radx_new_fname.replace('uf', 'nc')
		#radx_new_fname = 'SEA%s'%(radx_new_fname)

				# Here we need to run DROPS on the uf file created by RadxConvert
		drops_cmd = 'DROPS_v5 %s/realtime_drops/%s/%s -c %s -m VHS -o %s/realtime_drops/%s'%(base_path, stripped_radxdate, radx_old_fname, 
						drops_config_file, base_path, radx_new_fname)

		# now run the DROPS with the new and conforming filename
		os.system(drops_cmd)
		#print 'RAN DROPS'

		# now there's a new uf file with DROPS stuff, now need to convert that back into a netcdf with RadxConvert
		# First, get the most recent file in the directory

		drops_uf_file = os.path.basename(gtools.recent_file_glob('%s/realtime_drops'%(base_path), before_wildcard='*%s'%(stripped_month)))
		#print 'drops uf file: {}'.format(drops_uf_file)

		# Now use that filename to run radxconvert again in the same directory
		radx_conv_cmd = 'RadxConvert -f %s/realtime_drops/%s -outdir %s/realtime_drops'%(base_path, drops_uf_file, base_path)

		os.system(radx_conv_cmd)

		rx_path = '%s/realtime_drops/%s'%(base_path, stripped_radxdate)
		#print 'Converted DROPS output file'

		pd_old_fname = os.path.basename(gtools.recent_file_glob(rx_path, before_wildcard='*%s'%(stripped_month)))
		#print 'Post drops old filename: {}'.format(pd_old_fname)

		# okay, then the output file by the RadxConvert is a wonky filename, so move that output file to where it needs to be
		# gonna move it from 
		post_drops_prefix = 'cfrad.'

		pd_prefix_ind = pd_old_fname.index(post_drops_prefix) + len(post_drops_prefix)
		pd_radx_fbase = pd_old_fname[pd_prefix_ind:-1*(len(file_type)+1)]


		pd_date = pd_radx_fbase[:8]
		# pd_new_fname = pd_old_fname.replace(pd_radx_fbase, stripped_fbase)
		# pd_new_fname = pd_new_fname.replace('uf', 'nc')
		# pd_new_fname = 'SEA%s'%(pd_new_fname)

		pd_new_fname = 'SEA%s'%(radx_new_fname).replace('uf', 'nc')

		move_cmd = 'mv %s/%s %s/realtime_drops/%s'%(rx_path, pd_old_fname, base_path, pd_new_fname)
		os.system(move_cmd)


	else:
		print 'file: %s has already been gridded'%(stripped_fb)











