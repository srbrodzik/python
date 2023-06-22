# Contents
# assign_conv_radius_km
# radial_distance_mask
# background_intensity
# conv_core_cos_scheme
# make_initial_convsf_map
# incorporate_conv_radius
# write_netcdf_file_xyzt
# add_to_netcdf_file

def assign_conv_radius_km(bkgrnd_val,max_conv_radius,dbz_for_max_conv_radius):

    import sys

    # PURPOSE: from SHY algorithm based on the reflectivity of the background 
    # assigns the convective radius in KILOMETERS
    # for SHY algorithm as published
    # max_conv_radius = 5 km i.e convcore pixel plus 2 pixels if 2 km size pixels
    # dbz_for_max_conv_radius = 40
    # note as coded, which is similar to Steiner code,
    # pretty much assumes that max_conv_radius = 5 km

    # INPUTS:
    #       bkgrnd_val              = average refl val of background pixels
    #       max_conv_radius         = maximum convective radius (in km)
    #       dbz_for_max_conv_radius = min core dBZ for max conv radius mask
    # OUTPUTS:
    #       conv_radius_km          = convective radius in km

    if bkgrnd_val >= dbz_for_max_conv_radius:
        conv_radius_km=maxConvRadius
    elif bkgrnd_val >= dbz_for_max_conv_radius-5:
        conv_radius_km=maxConvRadius-1
    elif bkgrnd_val >= dbz_for_max_conv_radius-10:
        conv_radius_km=maxConvRadius-2
    elif bkgrnd_val >= dbz_for_max_conv_radius-15:
        conv_radius_km=maxConvRadius-3
    else:
        conv_radius_km=maxConvRadius-4

    return conv_radius_km

def radial_distance_mask(min_radius,max_radius,x_dim,y_dim,x_res,y_res):

    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    # PURPOSE:
    #       Creates either a circular mask for cartesian grid

    # INPUTS:
    #       min_radius = min radius of mask in km
    #       max_radius = max radius of mask in km
    #       x_dim      = size of x dimension
    #       y_dim      = size of y dimension
    #       x_res      = spacing between grid points in x direction
    #       y_res      = spacing between grid points in y direction

    # OUTPUTS:
    #       mask      = 2d cartesian mask

    debug = 1
    
    x_start = -(x_dim/2)+(x_res/2)   
    x_end   = (x_dim/2)-(x_res/2)   
    y_start = -(y_dim/2)+(y_res/2)   
    y_end   = (y_dim/2)-(y_res/2)

    x0 = x_start + x_end
    y0 = y_start + y_end
    
    x = np.linspace(x_start, x_end, x_dim)
    y = np.linspace(y_start, y_end, y_dim)

    x, y = np.meshgrid(x, y)
    r = np.sqrt((x - x0)**2 + (y - y0)**2)

    mask = np.logical_and(r >= min_radius, r < max_radius)
    #temp = np.ones(x.shape)
    out_array = np.ma.masked_where(~mask,temp)

    if debug:
        print >>sys.stderr, out_array

        fig, ax = plt.subplots()
        ax.set(xlabel='X', ylabel='Y', aspect=1.0)

        ax.scatter(x[mask], y[mask])

        plt.show()

    return out_array

def background_intensity(refl_array,x_dim,y_dim,win_diam,
                         kernel,missing_val):

    import numpy as np
    import math
    from scipy.signal import convolve2d as convolve

    # INPUTS:
    #       in subrout    in main                definition
    #       refl_array  = refl_2d_array        = 2d reflectivity
    #       x_dim       = x_dim                = x dimension of refl_array
    #       y_dim       = y_dim                = y dimension of refl_array
    #       win_diam    = max_bdg_win_diameter = num pixels in mask diam
    #       kernel      = win_mask_array       = background mask
    #       missing_val = missing_val          = refl missing value
    # OUTPUTS:
    #       bg_array    = background averages

    # convert non-missing values in refl_array to linear values
    # refl_array is masked so no need to hunt for good indices
    temp = np.ones(refl_array.shape)
    good_vals = np.ma.masked_where(refl_array.mask,temp)
    lin_array = pow(10.,refl_array/10)
    lin_sums = convolve(lin_array,kernel,mode='same',boundary='fill')
    num_pixs = convolve(good_vals,kernel,mode='same',boundary='fill')
    lin_avg = lin_sums/num_pixs
    bg_array = 10.0 * np.log10(lin_avg)
    
    return bg_array

def add_to_netcdf_file(ncid,convsf_array,missing_val):

    convsf_id = ncid.createVariable('convsf',nc.float,('time','y0','x0'),
                                    zlib=True,fill_value=missing_val)
    convsf_id.units = 'none'
    convsf_id.long_name = "Rain Type"
    convsf_id.stratiform = 1
    convsf_id.convective = 2
    convsf_id.other = 3
    convsf_id[:,:,:] = convsf_array

    ncid.close()
