# PURPOSE:
#  Determine the type of echo (rain_type) at each pixel of the gridded
#  radar data pattern at a chosen level:
#           - no surface rainfall  -> no-rain    = 0
#           - surface rainfall     -> stratiform = 1
#                                  -> convective = 2
#           - weak echo                            3

import sys
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import netCDF4 as nc4
import math

from conv_stra_routines import radial_distance_mask
from conv_stra_routines import background_intensity

input_dir = '/home/storm/brodzik/python/convsf/uw/test_refl_files'
dates = ['20111016']
in_prefix = 'radar'
field_name = 'REFL'
level = 5
output_dir = '/home/storm/brodzik/python/convsf/uw/test_convsf_files'
out_prefix = 'convsf_'
abs_conv_thresh = 64
min_Z_diff = 10
trunc_Z_conv_thresh = 40
dbz_for_max_conv_radius = 40
apply_rad_mask = 1
min_radius = 12
min_dbz_use = -50
weak_echo_thresh = 5
bkgrnd_radius = 11
max_conv_radius = 5
bkgrnd_pixel_flag = 1
debug = 1

for date in dates:
    for file in os.listdir(input_dir+'/'+date):
        if file.startswith(in_prefix):
            if debug:
                print >> sys.stderr, "file = ", file

            # Read file
            ncid = nc4.Dataset(input_dir+'/'+date+'/'+file,'a')
            time_array = ncid.variables['time'][:]
            x0_array = ncid.variables['x0'][:]
            y0_array = ncid.variables['y0'][:]
            z0_array = ncid.variables['z0'][:]
            refl_array = ncid.variables[field_name][:]
            missing_val =  ncid.variables[field_name]._FillValue
            #ncid.close()

            # Get file dimensions & horizontal resolution
            (t_dim,z_dim,y_dim,x_dim) = refl_array.shape
            x_res = x0_array[1] - x0_array[0]
            y_res = y0_array[1] - y0_array[0]
            max_radius = (x_dim/2)*x_res

            # Get level of refl to use
            refl_2d_array = np.squeeze(refl_array[:,level,:,:])

            # Declare one-D arrays
            # border_mask_array = np.array([x_dim,y_dim])
            # convsf_array = np.zeros([x_dim,y_dim]) + missing_val
            # bkgrnd_array = np.zeros([x_dim,y_dim]) + missing_val
            # conv_core_array = np.zeros([x_dim,y_dim])
            
            # Create background window mask
            # Initialize
            #   bkgrnd_radius in km; determine size of win_mask in pixels
            #   bdg_win_diameter must be odd to have pixel centered
            max_bdg_win_diameter = int( math.floor((float(bkgrnd_radius)/x_res)*2) )
            if math.floor(max_bdg_win_diameter/2) == max_bdg_win_diameter/2.0:
                max_bdg_win_diameter = max_bdg_win_diameter + 1
            win_mask_array = radial_distance_mask(0,bkgrnd_radius,max_bdg_win_diameter,
                                                  max_bdg_win_diameter,x_res,y_res)

            # CONVECTIVE MASKS
            
            # Initialize
            #   max_conv_radius in km; figure out max size in pixels
            #     of convective regions around conv core
            #   max_conv_diameter must be odd to have pixel centered
            max_conv_diameter = int( math.floor((float(max_conv_radius)/x_res)*2) )
            if math.floor(max_conv_diameter/2) == max_conv_diameter/2.0:
                max_conv_diameter = max_conv_diameter + 1
            conv_mask_array_all = ma.masked_array(np.ones([max_conv_radius,max_conv_diameter,max_conv_diameter]),mask
            mx = ma.masked_array(x, mask=[0, 0, 0, 1, 0])
            # Create (NOTE: 2nd param is radius in KM)
            for iradius in range(0,max_conv_radius):
                tmp_conv_mask = radial_distance_mask(0,iradius+1,max_conv_diameter,max_conv_diameter,x_res,y_res)
                conv_mask_array_all[iradius,:,:] = tmp_conv_mask
                print >>sys.stderr, "conv_mask_array_all[",iradius,",:,:] = ",conv_mask_array_all[iradius,:,:]
            if debug:
                print >>sys.stderr, "Done creating convective masks for all radii"

            # Compute background intensities
            bkgrnd_array = background_intensity(refl_2d_array,x_dim,y_dim,max_bdg_win_diameter,
                                                win_mask_array,missing_val)
            if debug:
                print >>sys.stderr, "Done calculating background intensities"

            # Identify conv cores
            conv_core_array = conv_core_cos_scheme(refl_2d_array,bkgrnd_array,x_dim,y_dim,
                                                   min_Z_diff,abs_conv_thresh,trunc_Z_conv_thresh,
                                                   bkgrnd_pixel_flag,missing_val,CS_CORE)
            if debug:
                print >>sys.stderr, "Done identifying convective cores"
            
            # Do initial assignment of convsf map & count conv cores
            make_initial_convsf_map(refl_2d_array,conv_core_array,convsf_array,
                                    weak_echo_thres,min_dbz_use,CS_CORE,CONV,
                                    WEAK_ECHO,NO_SFC_ECHO,SF,core_count,missing_val)
            if debug:
                print >>sys.stderr, "Done with initial assignment of convsf map"

            # REDO THIS SECTION USING ARRAY OPERATIONS
            # Go through refl_2d_array second time to apply conv radius to each core
            conv_rad_cnt = make_array(max_conv_radius,/integer,value=0)
            for j in range(0,y_dim):
                for i in range(0,x_dim):
                    index = (x_dim*j)+i             
                    if conv_core_array[index] == CS_CORE:
                        # Get int value of convRadius so can use as array index
                        conv_radius_km = assign_conv_radius_km(bkgrnd_array[index],
                                                               max_conv_radius,
                                                               dbz_for_max_conv_radius)
                        conv_radius = fix(floor(conv_radius_km))
                        # Define conv_mask_array and keep track of how many with each radius
                        if conv_radius >= 1 and conv_radius <= max_conv_radius:
                            conv_mask_array = conv_mask_array_all[*,*,convRadius-1]
                            conv_rad_cnt[conv_radius-1]=conv_rad_cnt[conv_radius-1]+1
                        elif:
                            print >>sys.stderr, "Unanticipated conv_mask_array size, conv_radius=" ,conv_radius
                            print >>sys.stderr, ". . . Need to change code for rad > ", max_conv_radius, ". Exiting"
                            return
                        # Modify convsf_array to incorporate points within radius
                        incorporate_conv_radius(convsf_array,i,j,x_dim,y_dim,
                                                max_conv_diameter,conv_mask_array)
            if debug:
                print >>sys.stderr, "Done going through array second time to assign conv radii to cores"

            # Create border mask
            # If apply_rad_mask set, cut out the corners
            # Else put in a border of missing_vals on the map since the algorithm
            #   cannot calculate bkgrnd Z properly at edges of map
            #   Use in cases where data beyond region of map
            if apply_rad_mask:
                max_map_radius = (x_dim * x_dim)/2.0
            elif:
                max_map_radius = (x_dim * x_res)/2.0 - bkgrnd_radius
            if debug:
                print >>sys.stderr, "Applying radial mask out to range = ", max_map_radius
            radial_distance_mask(min_radius,max_map_radius,x_dim,y_dim,
                                x_res,y_res,1,0,border_mask_array)

            # Apply border mask
            indMask = where(borderMaskArray ne 1,countMask)
            convsfArray[indMask] = missingVal
            if debug:
                print >>sys.stderr, "Done creating and applying border mask"

            # Output netcdf file
            # get output file name
            output_fn = output_dir+'/'+out_prefix+file
        
            add_to_netcdf_file(ncid,convsf_array,missing_val)
            if debug:
                print >>sys.stderr, "Done adding convsf array to netcdf file"
           
#---------------------------------------------------------------------------------
#x_dim = 400
#y_dim = 400
#x_res = 1.0
#y_res = 1.0
#max_radius = x_dim/2

# Test conv_radius_km routine
#conv_radius = assign_conv_radius_km(bkgrndVal,maxConvRadius,dBZforMaxConvRadius)
#print >>sys.stderr, "conv_radius = ", conv_radius

# Test radial_distance_mask routine
#mask = radial_distance_mask(min_radius,max_radius,x_dim,y_dim,x_res,y_res)
#print >>sys.stderr, "mask.shape = ", mask.shape
#print >>sys.stderr, "mask[0,0] = ", mask[0,0]
#print >>sys.stderr, "mask[100,100] = ", mask[100,100]
#print >>sys.stderr, "mask[200,200] = ", mask[200,200]
 
# Test background_intensity routine
#bg_array = background_intensity(refl_array,x_dim,y_dim,win_diam,mask)
#print >>sys.stderr, "bg_array.shape = ", bg_array.shape
#---------------------------------------------------------------------------------

