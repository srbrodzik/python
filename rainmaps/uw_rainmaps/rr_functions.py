import numpy as np
import math

def assign_conv_stra_mixed(rt_in, mixed_type, missing_value):

    """
    Description: Converts raintypes 1-6 from Powell et al algorithm to raintypes
    ranging from 1 to 3.  The mixed category (raintype=3) is assigned the value 
    defined by the input parameter 'mixed_type'.

    Inputs:
    rt_in - 2D array of raintypes; values range from 1 to 6 plus a missing val
       0 = no echo
       1 = stratiform
       2 = convective
       3 = mixed
       4 = isolated convective core (ICC)
       5 = isolated convective fringe (ICF)
       6 = weak echo
    missing_value - used in rt_in

    Outputs:
    rt_out - 2D array of raintype distilled to 3 types (stra, conv, mixed)
       1 = stratiform (includes rt_in values of 1, 5, and sometimes 3)
       2 = convective (includes rt_in values of 2, 4 and sometimes 3)
       3 = mixed (includes rt_in values of 6 and sometimes 3)
    rt_mask - locations of missing values

    """

    # Create new raintype array
    rt_out = np.zeros(np.shape(rt_in))
    rt_out[...] = -9999
    rt_out[rt_in==0] = 0
    rt_out[rt_in==1] = 1
    rt_out[rt_in==2] = 2
    rt_out[rt_in==3] = mixed_type
    rt_out[rt_in==4] = 2
    rt_out[rt_in==5] = 1
    rt_out[rt_in==6] = 0

    # Assign values to mask array
    rt_mask = np.zeros(np.shape(rt_out),dtype=np.int)
    rt_mask[(rt_out==missing_value)] = 1

    # Assign NaN to missing values of rt_out
    #rt_out[rt_mask==1] = np.nan
    
    return rt_out, rt_mask

def radial_distance_mask(min_radius, max_radius, xdim, ydim, x_spacing, y_spacing):

    """
    Description: Creates a radar coverage mask for a square grid

    Inputs:

    Outputs:

    """
    mask = np.zeros(shape=(ydim,xdim))
    
    center_pixel_y = int(np.floor(ydim/2))
    center_pixel_x = int(np.floor(xdim/2))

    for j in range(0,ydim):
        for i in range(0,xdim):
            y_range_sq = math.pow( ((center_pixel_y-j)*y_spacing), 2 )
            x_range_sq = math.pow( ((center_pixel_x-i)*x_spacing), 2 )
            dist = math.sqrt(x_range_sq+y_range_sq)
            if dist <= max_radius and dist >= min_radius:
                mask[j,i] = 1

    return mask
