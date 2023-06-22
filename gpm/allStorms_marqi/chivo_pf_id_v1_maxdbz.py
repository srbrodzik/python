# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 10:57:12 2021
Updated on Mon Mar 7 14:10:52 2022

@author: mrocq

PF code which fits an ellipse around a PF and then identifies TRMM features within the larger ellipse
Based on code from Weixin Xu, Kyle Chudler, Brody Fuchs, Zach Bruick, and Marqi Rocque

Current outputs are PF location (center of ellipse), area, conv/strat areas, mean dBZ
at each level, max dBZ at each level, and echo top height.


For 1km gridded CHIVO data from RELAMPAGO
"""

#############
import numpy as np  
from scipy.ndimage import label, generate_binary_structure, center_of_mass
import matplotlib.patches as patches
from numpy.linalg import eig
from matplotlib.lines import Line2D
#############

#takes in x and y points and returns parameters which define a ellipse that surrounds these points
#(Medioni et al. 2000; Nesbitt et al. 2006)
def get_ellipse(xs,ys):
    npts = len(xs)
    i11 = np.sum(ys**2) / npts
    i22 = np.sum(xs**2) / npts
    i12 = -np.sum(xs * ys) / npts
    tensor = [ [i11, i12], [i12,i22] ]
    eig_vals, eig_vecs = eig(tensor)
    semimajor = np.sqrt(eig_vals[0]) * 2.0
    semiminor = np.sqrt(eig_vals[1]) * 2.0
    major = semimajor * 2.0
    minor = semiminor * 2.0
    eig_vecs = eig_vecs[:,0]
    orientation = np.math.atan2(eig_vecs[1], eig_vecs[0]) - np.deg2rad(90.0)
    
    return major, minor, orientation


def PF_finder_max(refl, x, y, z, cs, hid, iwp, work_alt = 2.5, thresh_dbz = 0, dbz_conv = 40, dcc_hgt = 10, min_conv_n = 1000, min_strat_n = 40000):
  
    #set all nans to 0
    refl = np.ma.filled(refl,0)
    
    #set all negative values to 0
    refl[refl < 0] = 0
    
    #get horizontal grid spacing
    dx = np.diff(x)[0]

    #index of desired working altitude
    z_ind = np.where(z == work_alt)[0][0]
    
    #pull out reflecitvity at desired altitude
    refl_z = refl[z_ind]
    
    #turn x/y/z into 3D arrays
    y1,z1,x1 = np.meshgrid(y,z,x)
    
    #turn x/y into 2D arrays
    x,y = np.meshgrid(x,y)

    #get convective/stratiform components
    conv = np.zeros([x.shape[0], x.shape[1]])
    conv[np.where(cs > 0)] = 1
    conv[np.where(cs == 0)] = 0
    
#     conv[np.where(cs == 38)] = 1
#     conv[np.where(cs == 36)] = 1
#     conv[np.where(cs == 34)] = 1
#     conv[np.where(cs == 32)] = 1
#     conv[np.where(cs == 25)] = 1 #also include mixed in conv. class
    
    conv_inds = conv > 0
    
    strat = np.zeros([x.shape[0], x.shape[1]])
    strat[np.where(cs == 0)] = 1
    strat[np.where(cs > 0)] = 0
    
#     strat[np.where(cs == 14)] = 1
#     strat[np.where(cs == 16)] = 1
#     strat[np.where(cs == 18)] = 1
    
    strat_inds = strat > 0
                
    #get the composite reflectivity between 3 and 5 km (5 = 3km, 9 = 5km)
    max_dbz = np.zeros([x.shape[0],x.shape[1]])
    for i in range(y.shape[0]):
        for j in range(x.shape[0]):
            max_dbz[i,j] = np.nanmax(refl[5:10,i,j])


    #get the indices where max_dbz > dbz_conv and convective pixels exist
    conv_ind_40 = np.zeros([x.shape[0],x.shape[1]])
    
    for i in range(y.shape[0]):
        for j in range(x.shape[0]):
            if max_dbz[i,j] >= dbz_conv and conv_inds[i,j] == 1:
                conv_ind_40[i,j] = 1
            else:
                conv_ind_40[i,j] = 0
                
                
    #35 dBZ echo volume for lightning parameterization
    z35_ind = np.zeros([x1.shape[0],x1.shape[1],x1.shape[2]])
    
    for i in range(x1.shape[1]):
        for j in range(x1.shape[2]):
            for k in range(x1.shape[0]):
                if refl[k,i,j] >= 35:
                    z35_ind[k,i,j] = 1
                else:
                    z35_ind[k,i,j] = 0
                    
    #graupel volume for lightning parameterization    
    graup_ind = np.zeros([x1.shape[0],x1.shape[1],x1.shape[2]])
    
    for i in range(x1.shape[1]):
        for j in range(x1.shape[2]):
            for k in range(x1.shape[0]):
                if hid[k,i,j] == 7:    #graupel class is 12 in DROPS, is 7 and 8 in CSU
                    graup_ind[k,i,j] = 1
                elif hid[k,i,j] == 8:
                    graup_ind[k,i,j] = 1
                else:
                    graup_ind[k,i,j] = 0
                    
                    
    #hail index for hid percentage
    hail_ind = np.zeros([x1.shape[0],x1.shape[1],x1.shape[2]])
    
    for i in range(x1.shape[1]):
        for j in range(x1.shape[2]):
            for k in range(x1.shape[0]):
                if hid[k,i,j] == 9:    #hail class is 10,11 in DROPS, is 9 in CSU
                    hail_ind[k,i,j] = 1
#                 elif hid[k,i,j] == 11:
#                     hail_ind[k,i,j] = 1
                else:
                    hail_ind[k,i,j] = 0
                    
    #rain index for hid percentage
    liquid_ind = np.zeros([x1.shape[0],x1.shape[1],x1.shape[2]])
    
    for i in range(x1.shape[1]):
        for j in range(x1.shape[2]):
            for k in range(x1.shape[0]):
                if hid[k,i,j] == 1:    #rain indices 6-9 in DROPS, is 1, 2, 10 in CSU
                    liquid_ind[k,i,j] = 1
                elif hid[k,i,j] == 2:
                    liquid_ind[k,i,j] = 1
                elif hid[k,i,j] == 10:
                    liquid_ind[k,i,j] = 1
#                 elif hid[k,i,j] == 9:
#                     liquid_ind[k,i,j] = 1
                else:
                    liquid_ind[k,i,j] = 0
                    
    #snow index for hid percentage
    snow_ind = np.zeros([x1.shape[0],x1.shape[1],x1.shape[2]])
    
    for i in range(x1.shape[1]):
        for j in range(x1.shape[2]):
            for k in range(x1.shape[0]):
                if hid[k,i,j] == 4:    #snow indices 13,14 in DROPS, is 4,5 in CSU
                    snow_ind[k,i,j] = 1
                elif hid[k,i,j] == 5:
                    snow_ind[k,i,j] = 1
                else:
                    snow_ind[k,i,j] = 0
                    
                    
    #ice index for hid percentage
    ice_ind = np.zeros([x1.shape[0],x1.shape[1],x1.shape[2]])
    
    for i in range(x1.shape[1]):
        for j in range(x1.shape[2]):
            for k in range(x1.shape[0]):
                if hid[k,i,j] == 3:    #ice indices 15,16 in DROPS, is 3,6 in CSU
                    ice_ind[k,i,j] = 1
                elif hid[k,i,j] == 6:
                    ice_ind[k,i,j] = 1
                else:
                    ice_ind[k,i,j] = 0
        
    #masking array for finding groups. Points will be considered contiguous if they
    #touch a pixel on any side, including diagonally
    label_mask = generate_binary_structure(2,2)
    label_mask3 = generate_binary_structure(3,3)

    #group pixels that are adjacent
    pf_groups, n_groups = label(max_dbz > thresh_dbz, structure = label_mask)
    pf_groups_conv, n_groups_conv = label(conv_ind_40, structure = label_mask)
    pf_groups_strat, n_groups_strat = label(strat_inds, structure = label_mask)
    pf_groups_z35, n_groups_z35 = label(z35_ind, structure = label_mask3)
    pf_groups_graup, n_groups_graup = label(graup_ind, structure = label_mask3)

    #subtract 1 from pf_groups so that the indexing works out nicer
    pf_groups -= 1
    pf_groups_conv -= 1
    pf_groups_strat -= 1
    pf_groups_z35 -= 1
    pf_groups_graup -= 1

    #find the center of mass of the features
    pf_locs = center_of_mass(max_dbz >= thresh_dbz, pf_groups , np.arange(n_groups))
    pf_locs = [(x[int(l[0]),int(l[1])], y[int(l[0]),int(l[1])]) for l in pf_locs]
    
    pf_locs_conv = center_of_mass(conv_ind_40, pf_groups_conv , np.arange(n_groups_conv))
    pf_locs_conv = [(x[int(l[0]),int(l[1])], y[int(l[0]),int(l[1])]) for l in pf_locs_conv]
    
    pf_locs_strat = center_of_mass(strat_inds, pf_groups_strat , np.arange(n_groups_strat))
    pf_locs_strat = [(x[int(l[0]),int(l[1])], y[int(l[0]),int(l[1])]) for l in pf_locs_strat]    
    
    pf_locs_z35 = center_of_mass(z35_ind, pf_groups_z35, np.arange(n_groups_z35))
    pf_locs_z35 = [(z1[int(l[0]), int(l[1]), int(l[2])], y1[int(l[0]),int(l[1]),int(l[2])], x1[int(l[0]),int(l[1]),int(l[2])]) for l in pf_locs_z35] 

    pf_locs_graup = center_of_mass(graup_ind, pf_groups_graup, np.arange(n_groups_graup))
    pf_locs_graup = [(z1[int(l[0]), int(l[1]), int(l[2])], y1[int(l[0]),int(l[1]),int(l[2])], x1[int(l[0]),int(l[1]),int(l[2])]) for l in pf_locs_graup] 

    
    #set a prelimiary PF threshold with area > 4km2 and echotop height > 3km
    new_groups = []
    
    for group_num in np.arange(n_groups):
        
        pf_inds = pf_groups == group_num
            
        npix = np.sum(pf_inds)
        n_conv_pix  = np.sum((pf_inds) & (conv_inds))
        n_strat_pix  = np.sum((pf_inds) & (strat_inds))
        
        max_refl_by_alt_pf  = np.zeros(z.shape)
        for zi in range(len(z)):
            rain_inds = refl[zi, pf_inds] > 0
            if np.sum(rain_inds) > 0:
                max_refl_by_alt_pf[zi]  = np.max(refl[zi,pf_inds][rain_inds])

        echo_top_ind = len(max_refl_by_alt_pf) - np.argmax((max_refl_by_alt_pf > 0)[::-1]) - 1
        #echo_top.append(z[echo_top_ind])
    
        if (npix*dx**2) > 4 and z[echo_top_ind] > 3:
            new_groups.append(group_num)
            
    new_groups = np.array(new_groups)
   
    #define empty lists for params to be appended to
    ellipses           = []
    xloc               = []
    yloc               = []
    area               = []
    conv_area          = []
    strat_area         = []
    mean_refl_by_alt   = []
    max_refl_by_alt    = []
    echo_top           = []
    echo_top_40        = []
    mean_rr            = []
    max_rr             = []
    mean_conv_rr       = []
    mean_strat_rr      = []
    mean_rr_max        = []
    max_rr_max         = []
    mean_conv_rr_max   = []
    mean_strat_rr_max  = []
    per_graup          = []
    per_liqui          = []
    per_snow           = []
    per_ice            = []
    per_hail           = []
    iwp_max            = []
    iwp_mean           = []

    dwcc_ellipses = [[] for i in (new_groups)]
    dwcc_area = [[] for i in (new_groups)]
    dwcc_locx = [[] for i in (new_groups)]
    dwcc_locy = [[] for i in (new_groups)]
    dwcc_maxdbz_10km = [[] for i in (new_groups)]
    dwcc_mean_rr = [[] for i in (new_groups)]
    dwcc_mean_rr_max = [[] for i in (new_groups)]

    wcc_ellipses = [[] for i in (new_groups)]
    wcc_area = [[] for i in (new_groups)]
    wcc_locx = [[] for i in (new_groups)]
    wcc_locy = [[] for i in (new_groups)]
    wcc_mean_rr = [[] for i in (new_groups)]
    wcc_mean_rr_max = [[] for i in (new_groups)]
    
    dcc_ellipses = [[] for i in (new_groups)]
    dcc_area = [[] for i in (new_groups)]
    dcc_locx = [[] for i in (new_groups)]
    dcc_locy = [[] for i in (new_groups)]
    dcc_maxdbz_10km = [[] for i in (new_groups)]
    dcc_mean_rr = [[] for i in (new_groups)]
    dcc_mean_rr_max = [[] for i in (new_groups)]
    
    bsr_ellipses = [[] for i in (new_groups)]
    bsr_area = [[] for i in (new_groups)]
    bsr_locx = [[] for i in (new_groups)]
    bsr_locy = [[] for i in (new_groups)]
    bsr_mean_rr = [[] for i in (new_groups)]
    bsr_mean_rr_max = [[] for i in (new_groups)]
    
    z35_locx = [[] for i in (new_groups)]
    z35_locy = [[] for i in (new_groups)]
    z35_vol = [[] for i in (new_groups)]
    graup_locx = [[] for i in (new_groups)]
    graup_locy = [[] for i in (new_groups)]
    graup_vol = [[] for i in (new_groups)]
    graup_zmax = [[] for i in (new_groups)]

    i = 0
    for group_num in (new_groups):
        pf_inds = pf_groups == group_num
            
        npix = np.sum(pf_inds)
        n_conv_pix  = np.sum((pf_inds) & (conv_inds))
        n_strat_pix  = np.sum((pf_inds) & (strat_inds))
        
        xloc.append(pf_locs[group_num][0])
        yloc.append(pf_locs[group_num][1])

        area.append(npix*dx**2)
        conv_area.append(n_conv_pix*dx**2)
        strat_area.append(n_strat_pix*dx**2)


        pf_xs = x[pf_inds] - pf_locs[group_num][0]
        pf_ys = y[pf_inds] - pf_locs[group_num][1]

        major, minor, orientation = get_ellipse(pf_xs, pf_ys)

        ellipses.append(patches.Ellipse(pf_locs[group_num], width = major, height = minor, 
                                    angle = np.rad2deg(orientation), facecolor = 'None', edgecolor = 'k', lw = 2))


        mean_refl_by_alt_pf = np.zeros(z.shape)
        max_refl_by_alt_pf  = np.zeros(z.shape)
        graup_by_alt = np.zeros(z.shape)
        liqui_by_alt = np.zeros(z.shape)
        snow_by_alt = np.zeros(z.shape)
        ice_by_alt = np.zeros(z.shape)
        hail_by_alt = np.zeros(z.shape)
        all_dots = np.zeros(z.shape)
        
        for zi in range(len(z)):
            rain_inds = refl[zi, pf_inds] > 0
            if np.sum(rain_inds) > 0:
                mean_refl_by_alt_pf[zi] = np.mean(refl[zi,pf_inds][rain_inds])
                max_refl_by_alt_pf[zi]  = np.max(refl[zi,pf_inds][rain_inds])
                graup_by_alt[zi] = np.sum(graup_ind[zi,pf_inds][rain_inds])
                liqui_by_alt[zi] = np.sum(liquid_ind[zi,pf_inds][rain_inds])
                snow_by_alt[zi] = np.sum(snow_ind[zi,pf_inds][rain_inds])
                ice_by_alt[zi] = np.sum(ice_ind[zi,pf_inds][rain_inds])
                hail_by_alt[zi] = np.sum(hail_ind[zi,pf_inds][rain_inds])
                all_dots[zi] = np.sum(rain_inds)
                
        mean_refl_by_alt.append(mean_refl_by_alt_pf)
        max_refl_by_alt.append(max_refl_by_alt_pf)
        sum_dots = np.sum(all_dots)
        sum_graup = np.sum(graup_by_alt)
        sum_liquid = np.sum(liqui_by_alt)
        sum_snow = np.sum(snow_by_alt)
        sum_ice = np.sum(ice_by_alt)
        sum_hail = np.sum(hail_by_alt)
        
        per_graup.append(sum_graup/(sum_dots)*100)  
        per_liqui.append(sum_liquid/(sum_dots)*100)
        per_snow.append(sum_snow/(sum_dots)*100)  
        per_ice.append(sum_ice/(sum_dots)*100)
        per_hail.append(sum_hail/(sum_dots)*100)
        
        echo_top_ind = len(max_refl_by_alt_pf) - np.argmax((max_refl_by_alt_pf > 0)[::-1]) - 1
        echo_top.append(z[echo_top_ind])

        echo_top_40_ind = len(max_refl_by_alt_pf) - np.argmax((max_refl_by_alt_pf > dbz_conv)[::-1]) - 1
        
        if echo_top_40_ind > echo_top_ind:
            echo_top_40.append(0)
        else:
            echo_top_40.append(z[echo_top_40_ind])


        # add in rain rate, from DROPS2.0 paper (Chen et al. 2017, eq3):

        ze = 10**(refl_z/10) #using dbz at a specific altitude
        rr = 0.02*(ze**0.657)

        mean_rr.append(np.mean(rr[pf_inds]))
        max_rr.append(np.max(rr[pf_inds]))
        mean_conv_rr.append(np.mean(rr[(pf_inds) & (conv_inds)]))
        mean_strat_rr.append(np.mean(rr[(pf_inds) & (strat_inds)]))
        
        ze1 = 10**(max_dbz/10) #using max dbz in column
        rr1 = 0.02*(ze1**0.657)

        mean_rr_max.append(np.mean(rr1[pf_inds]))
        max_rr_max.append(np.max(rr1[pf_inds]))
        mean_conv_rr_max.append(np.mean(rr1[(pf_inds) & (conv_inds)]))
        mean_strat_rr_max.append(np.mean(rr1[(pf_inds) & (strat_inds)]))
        
        #get mean and max iwp of PF
        iwp_max.append(np.max(iwp[pf_inds]))
        iwp_mean.append(np.mean(iwp[pf_inds]))


        #look for TRMM features within the larger feature
        for conv_group_num in np.arange(n_groups_conv):

            if np.sum((pf_inds) & (pf_groups_conv == conv_group_num)) > min_conv_n/dx**2:

                inds = (pf_inds) & (pf_groups_conv == conv_group_num)
                z_10km = np.where(z == dcc_hgt)[0][0]
                tall_inds_10km = refl[z_10km, inds] > dbz_conv

                if np.sum(tall_inds_10km) > 0:
                    print ('DWCC found')

                    pf_xs = x[inds] - pf_locs_conv[conv_group_num][0]
                    pf_ys = y[inds] - pf_locs_conv[conv_group_num][1]

                    major, minor, orientation = get_ellipse(pf_xs, pf_ys)

                    dwcc_ellipses[i].append(patches.Ellipse(pf_locs_conv[conv_group_num],width = major, height = minor, 
                                      angle = np.rad2deg(orientation), facecolor = 'None', edgecolor = '#582f93', lw = 2))
                    dwcc_area[i].append(np.sum(inds)*dx**2)

                    dwcc_locx[i].append(pf_locs_conv[conv_group_num][0])
                    dwcc_locy[i].append(pf_locs_conv[conv_group_num][1])

                    dwcc_maxdbz_10km[i].append(np.max(refl[z_10km, inds]))
                    
                    ze = 10**(refl_z[inds]/10)
                    rr = 0.02*(ze**0.657)
                    dwcc_mean_rr[i].append(np.mean(rr))
                    
                    ze1 = 10**(max_dbz[inds]/10)
                    rr1 = 0.02*(ze1**0.657)
                    dwcc_mean_rr_max[i].append(np.mean(rr1))
    
                                           
                else:
                    print ('WCC found')

                    pf_xs = x[inds] - pf_locs_conv[conv_group_num][0]
                    pf_ys = y[inds] - pf_locs_conv[conv_group_num][1]

                    major, minor, orientation = get_ellipse(pf_xs, pf_ys)

                    wcc_ellipses[i].append(patches.Ellipse(pf_locs_conv[conv_group_num],width = major, height = minor, 
                                      angle = np.rad2deg(orientation), facecolor = 'None', edgecolor = 'dodgerblue', lw = 2))
                    wcc_area[i].append(np.sum(inds)*dx**2)
                    wcc_locx[i].append(pf_locs_conv[conv_group_num][0])
                    wcc_locy[i].append(pf_locs_conv[conv_group_num][1])
                    
                    ze = 10**(refl_z[inds]/10)
                    rr = 0.02*(ze**0.657)
                    wcc_mean_rr[i].append(np.mean(rr))
                    
                    ze1 = 10**(max_dbz[inds]/10)
                    rr1 = 0.02*(ze1**0.657)
                    wcc_mean_rr_max[i].append(np.mean(rr1))


            elif np.sum((pf_inds) & (pf_groups_conv == conv_group_num)) > 1/dx**2:

                inds = (pf_inds) & (pf_groups_conv == conv_group_num)

                z_10km = np.where(z == dcc_hgt)[0][0]
                tall_inds_10km = refl[z_10km, inds] > dbz_conv

                if np.sum(tall_inds_10km) > 0:
                    print ('DCC found')

                    pf_xs = x[inds] - pf_locs_conv[conv_group_num][0]
                    pf_ys = y[inds] - pf_locs_conv[conv_group_num][1]

                    major, minor, orientation = get_ellipse(pf_xs, pf_ys)

                    dcc_ellipses[i].append(patches.Ellipse(pf_locs_conv[conv_group_num],width = major, height = minor, 
                                      angle = np.rad2deg(orientation), facecolor = 'None', edgecolor = 'r', lw = 2))
                    dcc_area[i].append(np.sum(inds)*dx**2)
                    dcc_locx[i].append(pf_locs_conv[conv_group_num][0])
                    dcc_locy[i].append(pf_locs_conv[conv_group_num][1])

                    dcc_maxdbz_10km[i].append(np.max(refl[z_10km, inds]))
                    
                    ze = 10**(refl_z[inds]/10)
                    rr = 0.02*(ze**0.657)
                    dcc_mean_rr[i].append(np.mean(rr))
                    
                    ze1 = 10**(max_dbz[inds]/10)
                    rr1 = 0.02*(ze1**0.657)
                    dcc_mean_rr_max[i].append(np.mean(rr1))
                        
                                       
        #identify BSR       
        for strat_group_num in np.arange(n_groups_strat):
            if np.sum((pf_inds) & (pf_groups_strat == strat_group_num)) > min_strat_n/dx**2:
                print ('BSR found')

                inds = (pf_inds) & (pf_groups_strat == strat_group_num)

                pf_xs = x[inds] - pf_locs_strat[strat_group_num][0]
                pf_ys = y[inds] - pf_locs_strat[strat_group_num][1]

                major, minor, orientation = get_ellipse(pf_xs, pf_ys)

                bsr_ellipses[i].append(patches.Ellipse(pf_locs_strat[strat_group_num],width = major, height = minor, 
                                  angle = np.rad2deg(orientation), facecolor = 'None', edgecolor = 'orange', lw = 2))
                bsr_area[i].append(np.sum(inds)*dx**2)
                bsr_locx[i].append(pf_locs_strat[strat_group_num][0])
                bsr_locy[i].append(pf_locs_strat[strat_group_num][1])
                
                ze = 10**(refl_z[inds]/10)
                rr = 0.02*(ze**0.657)
                bsr_mean_rr[i].append(np.mean(rr))

                ze1 = 10**(max_dbz[inds]/10)
                rr1 = 0.02*(ze1**0.657)
                bsr_mean_rr_max[i].append(np.mean(rr1))


        #identify 35-dBZ volumes
        for z35_group_num in np.arange(n_groups_z35):
            inds = (pf_inds) & (pf_groups_z35 == z35_group_num)

            npix = np.sum(inds)
            
            if (npix*dx*dx*0.5) > 2:
                z35_vol[i].append(npix*dx*dx*0.5)
                z35_locx[i].append(pf_locs_z35[z35_group_num][2])
                z35_locy[i].append(pf_locs_z35[z35_group_num][1])


        #identify graupel volumes
        for graup_group_num in np.arange(n_groups_graup):
            inds = (pf_inds) & (pf_groups_graup == graup_group_num)

            npix = np.sum(inds)

            if (npix*dx*dx*0.5) > 2:
                graup_vol[i].append(npix*dx*dx*0.5)  
                zs = z1[inds]
                zmax = np.max(zs)
                graup_zmax[i].append(zmax)
                graup_locx[i].append(pf_locs_graup[graup_group_num][2])
                graup_locy[i].append(pf_locs_graup[graup_group_num][1])
                
        i = i+1
      
    #spit out all parameters
    return ellipses, xloc, yloc, area, conv_area, strat_area, mean_refl_by_alt, \
              max_refl_by_alt, echo_top, echo_top_40, mean_rr, max_rr, \
              mean_conv_rr, mean_strat_rr, mean_rr_max, max_rr_max, mean_conv_rr_max, \
              mean_strat_rr_max, iwp_max, iwp_mean, dcc_ellipses, dcc_area, \
              dcc_locx, dcc_locy, dcc_maxdbz_10km, dcc_mean_rr, dcc_mean_rr_max, \
              dwcc_ellipses, dwcc_area, dwcc_locx, dwcc_locy, dwcc_maxdbz_10km, \
              dwcc_mean_rr, dwcc_mean_rr_max, wcc_ellipses, wcc_area, \
              wcc_locx, wcc_locy, wcc_mean_rr, wcc_mean_rr_max, bsr_ellipses, bsr_area, \
              bsr_locx, bsr_locy, bsr_mean_rr, bsr_mean_rr_max, z35_locx, z35_locy, \
              z35_vol, graup_locx, graup_locy, graup_vol, graup_zmax, per_graup, \
              per_liqui, per_snow, per_ice, per_hail

# ##############################################################################