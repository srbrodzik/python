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

    mask = np.logical_and(r > min_radius, r < max_radius)

    if debug:
        fig, ax = plt.subplots()
        ax.set(xlabel='X', ylabel='Y', aspect=1.0)

        ax.scatter(x[mask], y[mask])

        plt.show()

    return mask

def background_intensity(refl_array,x_dim,y_dim,win_diam,
                         kernel,missing_val):

    import numpy as np
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

def conv_core_cos_scheme(refl_array,bg_array,x_dim,y_dim,
                         min_Z_diff,abs_conv_thresh,trunc_Z_conv_thresh,
                         bkgrnd_pixel_flag,missing_val,CS_CORE):

    # INPUTS:
    #       refl_array      = refl_2d_array   = 2d reflectivity
    #       bg_array        = bkgrnd_array    = background averages       
    #       x_dim           = x_dim           = x dimension of refl_array
    #       y_dim           = y_dim           = y dimension of refl_array
    #       min_Z_diff      = min_Z_diff      = value used in cosine function
    #       abs_conv_thresh = abs_conv_thresh = value used in cosine function
    #       trunc_Z_conv_thresh = trunc_Z_conv_thresh = refl threshold above
    #                                             which point is automatically
    #                                             a convective core
    #       bkgrnd_pixel_flag   = bkgrnd_pixel_flag   = 0=pixel, 1=bkgrnd
    #       missing_val     = missing_val     = refl missing value
    #       CS_CORE         = CS_CORE         = indicates conv core (set to 3)
    # OUTPUTS:
    #       conv_core_array = conv_core_array = convective cores
    
    HUGE = 1000.

    # if refl above threshold, then automatically a core
    core_ind = where(refl_array >= trunc_Z_conv_thresh,core_cnt)
    if core_cnt > 0:
        conv_core_array[core_ind] = CS_CORE

    # if refl below threshold, see if refl sufficiently above bkgrnd value
    z_diff = make_array(x_dim,y_dim,/float,val=HUGE)
    ind1 = where(refl_array != missing_val and refl_array < trunc_Z_conv_thresh and bg_array != missing_val and bg_array < 0,count1)
    if count1 > 0:
        z_diff[ind1] = min_Z_diff

    ind2 = where(refl_array != missing_val and refl_array < trunc_Z_conv_thresh and bg_array != missing_val and bg_array >= 0 and bkgrnd_pixel_flag == 1,count2)
    if count2 > 0:
        z_diff[ind2] = min_Z_diff * cos(!pi*bg_array[ind2]/(2*abs_conv_thresh)) 

    ind3 = where(refl_array != missing_val and refl_array < trunc_Z_conv_thresh and bg_array != missing_val and bg_array >= 0 and bkgrnd_pixel_flag == 0,count3)
    if count3 > 0:
        z_diff[ind3] = min_Z_diff * cos(!pi*refl_array[ind3]/(2*abs_conv_thresh))
  
    ind4 = where(refl_array != missing_val and refl_array < trunc_Z_conv_thresh and z_diff < 0,count4)
    if count4 > 0:
        z_diff[ind4] = 0
  
    ind5 = where(refl_array != missing_val and refl_array < trunc_Z_conv_thresh and bg_array != missing_val and (refl_array-bg_array) >= z_diff,count5)
    if count5 > 0:
        conv_core_array[ind5] = CS_CORE

    return conv_core_array

def make_initial_convsf_map(refl_array,conv_core_array,convsf_array,
                            weak_echo_thres,min_dbz_use,CS_CORE,CONV,
                            WEAK_ECHO,NO_SFC_ECHO,SF,core_count,missing_val):

    # INPUTS:
    #       refl_array       = refl_2d_array      = reflectivity array
    #       conv_core_array  = conv_core_array    = convective core array
    #       convsf_array     = convsf-array       = conv-stra array
    #       weak_echo_thresh = weak_echo_thresh   = weak echo threshold
    #       min_dbz_use      = min_dbz_use        = minimum dBZ for echo
    #       CS_CORE                               = value of conv core in convCoreArray (3)
    #       CONV                                  = value of convective pixel in convsfArray (2)
    #       WEAK_ECHO                             = value of weak echo pixel in convsfArray (3)
    #       NO_SFC_ECHO                           = value of no surface echo pixel in convsfArray (0)
    #       SF                                    = value of stratiform pixel in convsfArray (1)
    #       core_count       = core_count         = number of identified cores in convCoreArray
    #       missing_val      = missing_val        = reflectivity missing value
    # OUTPUTS:
    #       convsf_array     = convsf_array       = conv-stra array
    
    conv_core_ind = where(conv_core_array == CS_CORE,core_count)
    if coreCount > 0:
        convsf_array[conv_core_ind] = CONV

    no_echo_ind = where(refl_array != missing_val and refl_array < min_dbz_use,
                        noEchoCount)
    if noEchoCount > 0:
        convsf_array[noEchoInd] = NO_SFC_ECHO

    weak_echo_ind = where(refl_array != missing_val and refl_array >= min_dbz_use and refl_array < weak_echo_thresh,weak_echo_count)
    if weak_echo_count > 0:
        convsf_array[weak_echo_ind] = WEAK_ECHO
  
    stra_ind = where(refl_array != missing_val and convsf_array != CONV and convsf_array != WEAK_ECHO and convsf_array != NO_SFC_ECHO, stra_count)
    if stra_count > 0:
        convsf_array[stra_ind] = SF

    return convsf_array

def incorporate_conv_radius(convsf_array,x_index,y_index,x_dim,y_dim,max_conv_diam,mask):

    # INPUTS:
    #       convsf_array  = convsf array
    #       x_index       = x index in convsfArray
    #       y_index       = y index in convsfArray
    #       x_dim         = num pixels in convsfArray in x dimension
    #       y_dim         = num pixels in convsfArray in y dimension
    #       max_conv_diam = max conv diameter in pixels
    #       mask          = conv mask array for this pixel
    # OUTPUTS:
    #       convsf_array  = convsf array

    NO_SFC_ECHO = 0
    SF = 1
    CONV = 2
    WEAK_ECHO = 3

    # Initialize
    i_mask = 0
    j_mask = 0

    # NOTE: the conv_radius is always odd as it includes the conv_core 
    #       pixel at the center
    pixelRadius = fix( floor(max_conv_diam/2) )

    if x_index >= 0 and x_index < x_dim and y_index >= 0 and y_index < y_dim:
     
        # Loop through the points in the square window on the full-sized map

        if debug:
            print >>sys.stderr, "y_index-pixel_radius = ", y_index-pixel_radius
            print >>sys.stderr, "y_index+pixel_radius = ", y_index+pixel_radius

        for j in range(y_index-pixel_radius,y_index+pixelRadius+1):

            if debug:
                print sys.stderr, "x_index-pixel_radius = ", x_index-pixel_radius
                print sys.stderr, "x_index+pixel_radius = ", x_index+pixel_radius

            i_mask = 0
            
            if debug:
                print >>sys.stderr, "i_mask = ", i_mask, " and j_mask = ", j_mask

            for i in range(x_index-pixel_radius,x_index+pixel_radius+1):

                if debug:
                    print >>sys.stderr, "mask index = ", (max_conv_diam*j_mask)+i_mask
                    print >>sys.stderr, "convsf index = ", (x_dim*j)+i
                    print >>sys.stderr, "i = ", i, " and j = ", j
                    print >>sys.stderr, "x_dim = ", x_dim, " and y_dim = ", y_dim
                    print >>sys.stderr, "i_mask = ", i_mask, " and j_mask = ", j_mask
                    print >>sys.stderr, "max_conv_diam = ", max_conv_diam

                if i >= 0 and i < x_dim and j >= 0 and j < y_dim and mask[(max_conv_diam*j_mask)+i_mask] eq 1:

                    if convsf_array[(x_dim*j)+i] != NO_SFC_ECHO:
              
                        convsf_array[(x_dim*j)+i] = CONV
              
                i_mask = i_mask + 1

            j_mask = j_mask + 1

    return convsf_array

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
