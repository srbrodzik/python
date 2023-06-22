## ________________________________________________________________________________________________
'''

Convective & Stratiform Masking (separation, identification) algorithm for 2D gridded reflectivity from observations or models.

Name:
    conv_stra_mask

Purpose:
    This program takes a 2D gridded reflectivity field and designates each point convective or stratiform based upon either the [ Steiner, Houze and Yuter (SHY 1995 - JAMS) ] or [ Yuter and Houze (YH 1998 - QJMS; YH 1997 - JAM)  ] algorithm. The module can now handle either the Cartesian-gridded data (Lat-Lon coordinates), or the Lambert-conformal-gridded data (constant distance grid). The latter is common among models output (e.g. WRF CONUS runs).

Calling sequence:
    conv_stra_sep
    
Run: 
    CS_mask, Conv_Cent, Bkgnd_REFL = conv_stra_sep( dbz_data, XLAT_data, XLONG_data, grid_res, coor_type, method )

Input:
    
    dbz_data    2-D radar reflectivity data (dBZ)
    
    XLAT_data   2-D Latitude coordinates (e.g., XLAT from WRF model output)
    
    XLONG_data  2-D Longitude coordinates (e.g., XLONG from WRF model output)
                NOTE: The XLAT/XLONG are only used when the input data is Cartesian-gridded (in Lat-Lon coordinates)
    
    grid_res    Resolution of the grid (degree if data is Cartesian-gridded; km for Lambert-conformal-gridded data)
    
    coor_type   Types of data coordinates
                NOTE: 'C' - artesian-gridded data
                      'L' - Lambert-conformal-gridded data)
    
    method      Indicates whether to use the SHY or YH methodology for convective centers calculations.
                NOTE: 'SHY' - SHY (1995) used where radius around center following profile used in the paper (Fig. 7)
                      'YH'  - YH (1997, 1998) used where paramters can be tuned (Function B1 in YH97; 1 in YH98)

Optional Keywords:
    
    CoreThresh  Intensity threshold to identify convective centers. Reflectivity of at least 40 dBZ is automatically 
                labeled as a convective center (YH98 uses 46)
    
    tune_thresh Parameter background reflectivity threshold for Peakedness criteria in SHY95 (Function 2, Fig. 7)
    
    bg_diff     Parameter background reflectivity difference for Peakedness criteria in SHY95 (Function 2, Fig. 7)
    
    fill_dbz    Filling dBZ value for missing data points.
    
Output:
    
    CS_mask     2-D array with the Convective-Stratiform Masks produced by this algorithm based on Surrounding area 
                criteria (Fig. 6b in SHY95) according to the background reflectivity
                = NaN     missing data
                = 0       stratiform (where data is good & reflectivity >= 15 dBZ)
                = 1       convective (background reflectivity <= 25 dBZ)
                = 2       convective (25 < background reflectivity <= 30 dBZ)
                = 3       convective (30 < background reflectivity <= 35 dBZ)
                = 4       convective (35 < background reflectivity <= 40 dBZ)
                = 5       convective (background reflectivity > 40 dBZ)
    
    Conv_Cent   2-D array with the Convective Centers identified by this algorithm based on the Intensity and Peakedness
                criteria (SHY95) according to the actual and background reflectivity
                = -1      Not a Convective Center
                = 1       Convective Center identified by the Intensity criteria
                = 21      Convective Center identified by the Peakedness criteria 1 (Function 2 top in SHY95)
                = 22      Convective Center identified by the Peakedness criteria 2 (Function 2 middle in SHY95)
                = 23      Convective Center identified by the Peakedness criteria 3 (Function 2 bottom in SHY95)
    
    Bkgnd_REFL  2-D array of the Bcakground Reflectivity calculated from the actual reflectivity (dBZ).
                The calculation is done by a smoothing (moving averaging) filter. According to SHY95, a 11-km radius 
                circle was used to mimic the 400-km^2 area from Churchill and Houze (1984). Herein, we directly use the 
                20-km length (sqrt(400 km^2)) divided by the grid_res (4-km in this case) for the filter size.

Author and history:
    
    Apr 2011    Created by Nick Guy (CSU)
                Modified to speed up performance.
    Apr 2011    Corrected the mean reflectivity of background field
                to take mean of linear Z and not mean of dBZ.
    Dec 2011    Arrays initialized to Nan.
                Missing value added and updated in procedure.
    May 2018    Brody Fuchs
    Oct 2018    Modifed by Zach Bruick for use with WRF data for RELAMPAGO, utilizing Hopper and Schumacher (2012) 
                numbers for a, b, and min dBZ for convective classification.
    Dec 2021    Modified by Hungjui Yu for use with Hires WRF CONUS CTRL & PGW runs.
                Add functions to the algorithm that the module can now handle either the Cartesian-gridded data (Lat-
                Lon coordinates), or the Lambert-conformal-gridded data (constant distance grid). The latter is common 
                among models output (e.g. WRF CONUS runs).

References:
    Steiner, M., R. A. Houze, Jr., and S. E. Yuter, 1995: Climatological characterization of three-dimensional storm 
        structure from operational radar and rain gauge data. J. Appl. Meteor., 34, 1978-2007.
    Yuter, S. E., and R. A. Houze, Jr., 1997: Measurements of raindrop size distributions over the Pacific warm pool and 
        implications for Z-R relations. J. Appl. Meteor., 36, 847-867.
    Yuter, S. E., and R. A. Houze, Jr., 1998: The natural variability of precipitating clouds over the western Pacific 
        warm pool. Quart. J. Roy. Meteor. Soc., 124, 53-99.

'''
## ________________________________________________________________________________________________

import numpy as np
import scipy.ndimage as ndi
import os

## ________________________________________________________________________________________________

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    
    re = 6371.0
    
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    return c * re

## ________________________________________________________________________________________________

def pythagorean_dist(loc_id_x, loc_id_y, grid_length):
    """
    Calculate the distance between two points by their location id
    based on the Pythagorean Theorem on Euclidean (plain) geometry.
    """
    return np.sqrt( (loc_id_x*grid_length)**2 + (loc_id_y*grid_length)**2 )

## ________________________________________________________________________________________________

def mask_conv_area_C(bkgnd_range_id, outarray, latarray, lonarray, grid_res, dist_thresh, assign_value):
    
    ## Set the Box Size for processing speed:
    ## The box_size should be at least 5-km by the grid_res according to SHY95.
    box_size = np.ceil((5/(grid_res*110)))
    box_size = int(box_size)
    
    ## Loop through each Convective Center:
    for ii in range(0, len(bkgnd_range_id[0])):
        
        ## Get the Lat/Lon for the Convective Center:
        lat_ind = bkgnd_range_id[0][ii]
        lon_ind = bkgnd_range_id[1][ii]
        
        lat1 = latarray[lat_ind, lon_ind]
        lon1 = lonarray[lat_ind, lon_ind]

        ## Get just the box area (box_size) surrounding the Convective Center:
        i_lat_min = np.max([0, lat_ind-box_size])
        i_lat_max = np.min([latarray.shape[0], lat_ind+box_size])

        i_lon_min = np.max([0, lon_ind-box_size])
        i_lon_max = np.min([lonarray.shape[1], lon_ind+box_size])

        lat2 = latarray[i_lat_min:i_lat_max+1, i_lon_min:i_lon_max+1]
        lon2 = lonarray[i_lat_min:i_lat_max+1, i_lon_min:i_lon_max+1]

        ## Calculate the array of distance (in the box) to the convective center:
        d = haversine(lat1, lon1, lat2, lon2)

        ## Mask the surrounding Convective Area according to the distance threshold:
        close_d = np.where(d <= dist_thresh)
        
        if len(close_d[0]):
            close_d_global = ( close_d[0] + lat_ind - box_size
                             , close_d[1] + lon_ind - box_size
                             )
            outarray[close_d_global] = assign_value

    return outarray

## ________________________________________________________________________________________________

def mask_conv_area_L(bkgnd_range_id, outarray, grid_res, dist_thresh, assign_value):
    
    ## Initialize the location id arrays for Lat/Lon grids:
    locid_lat_0 = np.tile(np.arange(outarray.shape[0]),(outarray.shape[1],1))
    locid_lat_0 = np.transpose(locid_lat_0)
    locid_lon_0 = np.tile(np.arange(outarray.shape[1]),(outarray.shape[0],1))
    
    ## Set the Box Size for processing speed:
    ## The box_size should be at least 5-km by the grid_res according to SHY95.
    box_size = np.ceil(5/grid_res)
    box_size = int(box_size)
    
    ## Loop through each Convective Center:
    for ii in range(0, len(bkgnd_range_id[0])):
        
        ## Get the location id for the Convective Center:
        locid_lat = bkgnd_range_id[0][ii]
        locid_lon = bkgnd_range_id[1][ii]
        
        ## Adjust the location id arrays by setting the id for Convective Center as (0, 0):
        # locid_lat_adj = locid_lat_0 - locid_lat
        # locid_lon_adj = locid_lon_0 - locid_lon
        
        ## Get just the box area (box_size) surrounding the Convective Center:
        i_locid_lat_min = np.max([0, locid_lat-box_size])
        i_locid_lat_max = np.min([locid_lat_0.shape[0], locid_lat+box_size])

        i_locid_lon_min = np.max([0, locid_lon-box_size])
        i_locid_lon_max = np.min([locid_lon_0.shape[1], locid_lon+box_size])
        
        # locid_lat_box = locid_lat_adj[i_locid_lat_min:i_locid_lat_max+1, i_locid_lon_min:i_locid_lon_max+1] 
        # locid_lon_box = locid_lon_adj[i_locid_lat_min:i_locid_lat_max+1, i_locid_lon_min:i_locid_lon_max+1] 
        
        locid_lat_box_2 = locid_lat_0[i_locid_lat_min:i_locid_lat_max+1, i_locid_lon_min:i_locid_lon_max+1] - locid_lat
        locid_lon_box_2 = locid_lon_0[i_locid_lat_min:i_locid_lat_max+1, i_locid_lon_min:i_locid_lon_max+1] - locid_lon
        
        ## Calculate the array of distance (in the box) to the convective center:
        # dist_arr = np.sqrt( (locid_lat_adj*grid_res)**2 + (locid_lon_adj*grid_res)**2 )
        # dist_arr = np.sqrt( (locid_lat_box*grid_res)**2 + (locid_lon_box*grid_res)**2 )
        
        dist_arr = pythagorean_dist(locid_lat_box_2, locid_lon_box_2, grid_res)

        ## Mask the surrounding Convective Area according to the distance threshold:
        close_dist = np.where(dist_arr <= dist_thresh)
        
        if len(close_dist[0]):
            close_dist_global = ( close_dist[0] + locid_lat - box_size
                                , close_dist[1] + locid_lon - box_size
                                )
            outarray[close_dist_global] = assign_value

    return outarray
    
## ________________________________________________________________________________________________

def conv_stra_sep(dbz, lat2d, lon2d, grid_res, coor_type, method, CoreThresh=40.0, tune_thresh=42.43, bg_diff=10, fill_dbz=0.0):
    
    ## Function arguments settings:
    # coor_type = 'C' if coor_type is None else coor_type 
    # method = 'SHY' if method is None else method
    assert coor_type in ['C', 'L'] # make sure the coordinates type is either "C"artesian or "L"ambert-conformal.
    assert method in ['SHY', 'YH'] # make sure the method is either SHY or YH
    
    ## Initialize the mask fields:
    cs = np.full_like(dbz, np.nan) # Conv./Stra. Mask.
    # cc = np.zeros_like(cs) - 1     # Convective Centers.
    
    ## Set good/bad data booleans:
    ## The bad indicates missing dBZ data points or dBZ <= 0 
    ## (The smallest dBZ in models runs is -30.0; and the dBZ is simulated as S-band radar(?)):
    boo_bad = ((dbz == np.nan) | (dbz <= 0.0))
    boo_good = np.logical_not(boo_bad)
    
    ## Set the Conv./Stra. Masks anywhere data is good and dBZ greater than 15 dBZ (YH97, 98), give it a 0.
    ## 0 is stratiform:
    cs[((boo_good) & (dbz >= 0.0))] = 0
    
    ## Fill the missing grids with specified dBZ values:
    dbz[boo_bad] = fill_dbz

    ## Calculate a smoothed (averaged) background reflectivity field using a filter:
    ## According to SHY95, a 11-km radius circle was used to mimic the 400-km^2 area from Churchill and Houze (1984).
    ## Herein, we directly use the 20-km length (sqrt(400)) divided by the grid_res (4-km in this case) for the filter size.
            
    ## Calculate linear reflectivity:
    zlin = 10.**(dbz/10.)
    
    # bkgnd_lin = ndi.median_filter(zlin, size=, mode='nearest')
    bkgnd_lin = ndi.uniform_filter(zlin, size=20/grid_res, mode='nearest')
    bkgnd = 10.*np.log10(bkgnd_lin)
    bkgnd[boo_bad] = fill_dbz
    # bkgnd[(bkgnd == np.nan))] = fill_dbz
    
    ## Convective Centers identification 1: Intensity (SHY95):
    cc = np.where( (dbz >= CoreThresh), 1, -1 )
    
    ## Convective Centers identification 2: Peakedness (SHY95):
    ## (a): Function (2) top in SHY95.
    
    CCI2a = (bkgnd < 0.) & ((dbz - bkgnd) >= bg_diff)
    cc[CCI2a] = 21

    ## Convective Centers identification 2: Peakedness (SHY95):
    ## (b): Function (2) middle in SHY95.
    
    if method == 'SHY':
        CCI2b = (bkgnd >= 0.) & (bkgnd < tune_thresh) & (dbz-bkgnd >= (bg_diff-(bkgnd**2.)/180.))
    elif method == 'YH':
        # This line uses YH (1998) climatological tuning algorithm (a=10, b=100):
        CCI2b = (bkgnd >= 0.) & (bkgnd < tune_thresh) & (dbz-bkgnd > 10.*np.cos((np.pi*bkgnd)/(2.*100.)))
    else:
        raise ValueError('The method is not defined!')

    cc[CCI2b] = 22
    
    ## Convective Centers identification 2: Peakedness (SHY95):
    ## (c): Function (2) bottom in SHY95.

    CCI2c = (bkgnd >= tune_thresh)
    cc[CCI2c] = 23

    ## Convective Centers identification 3: Surrounding area (SHY95):
    ## Use the medium relation as in Fig.6b of SHY95.
    
    if coor_type == 'C':
    
        bkgnd_range_id = np.where((cc > 0) & (bkgnd <= 25.)) # & (bkgnd >= 15.) 
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_C(bkgnd_range_id, cs, lat2d, lon2d, grid_res, 1.0, 1)
            # cs[bkgnd_range_id] = 1

        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 25.) & (bkgnd <= 30.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_C(bkgnd_range_id, cs, lat2d, lon2d, grid_res, 2.0, 2)
            # cs[bkgnd_range_id] = 2

        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 30.) & (bkgnd <= 35.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_C(bkgnd_range_id, cs, lat2d, lon2d, grid_res, 3.0, 3)
            # cs[bkgnd_range_id] = 3

        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 35.) & (bkgnd <= 40.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_C(bkgnd_range_id, cs, lat2d, lon2d, grid_res, 4.0, 4)
            # cs[bkgnd_range_id] = 4

        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 40.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_C(bkgnd_range_id, cs, lat2d, lon2d, grid_res, 5.0, 5)
            # cs[bkgnd_range_id] = 5
    
    elif coor_type == 'L':
        
        bkgnd_range_id = np.where((cc > 0) & (bkgnd <= 25.)) # & (bkgnd >= 15.)  
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_L(bkgnd_range_id, cs, grid_res, 1.0, 1)

        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 25.) & (bkgnd <= 30.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_L(bkgnd_range_id, cs, grid_res, 2.0, 2)
            
        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 30.) & (bkgnd <= 35.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_L(bkgnd_range_id, cs, grid_res, 3.0, 3)
            
        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 35.) & (bkgnd <= 40.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_L(bkgnd_range_id, cs, grid_res, 4.0, 4)
            
        bkgnd_range_id = np.where((cc > 0) & (bkgnd > 40.))
        if len(bkgnd_range_id[0]):
            cs = mask_conv_area_L(bkgnd_range_id, cs, grid_res, 5.0, 5)
                        
    else:
        raise ValueError('The coordinates type is not defined!')
        
        
    return cs, cc, bkgnd

## ________________________________________________________________________________________________
## ________________________________________________________________________________________________
