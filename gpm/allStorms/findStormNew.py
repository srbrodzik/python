#!/usr/bin/python3

# Find contiguous pixels with data in the 3D reflectivity matrix 
# matching the thereshold for reflectivty . Identify individual storms!!

# INPUTS:
# refl_3D - 3d grid of reflectivity values

# OUTPUTS:
# id_storm - array containing sorted IDs of all unique storms
# npix_str - array containing pixel counts of each unique storm
# grid_storm - array of same size as refl_3D filled with storm IDs

import math
import numpy as np

# For testing
#refl_3D = np.ones((3,4,5))*-1
#refl_3D[0,2,2]=10
#refl_3D[0,3,2:5]=10
#refl_3D[1,1,2:5]=10
#refl_3D[1,2]=10
#refl_3D[0,0,0]=5

def findStorm(refl_3D):

    fillValue = -999

    #select only pixels with values greater than 0
    tam = refl_3D.shape
    cswath = np.ones(tam,dtype='int16') * fillValue

    # if there are values greater than zero in refl_3D, continue, else exit
    if len(cswath[refl_3D >0]) > 0:
        cswath[refl_3D > 0] = 0
    
        # expand the refl3D matrix 1 row,col,lev to either side
        grid = np.pad(cswath,(1,1),'constant',constant_values=(-999,-999))
        tam2 = grid.shape

        # objnum changes when a new storm is found
        objnum=1

        # this is evaluating from bottom to top (UpDown)
        for UD in range(0,tam2[0]-1):
            
            # this is evaluating from North to South  (lats-rows)
            for NS in range(0,tam2[1]-1):
                
                # this is evaluating from East to West (lons-cols)
                for EW in range(0,tam2[2]-1):

                    # print('UD=',UD,' NS=',NS,' EW=',EW,' grid[UD,NS,EW]=',grid[UD,NS,EW])
         
                    if grid[UD,NS,EW] != fillValue:
                        #grab the cells around the cell of interest in all directions
                        #NUMBERING OF SURROUNDING CELLS %%%
                        #    samelevel =  1   2   3    uplevel =  9  10  11    downlevel = 18  19  20    %
                        #                 4       5              12  13  14                21  22  23    %
                        #                 6   7   8              15  16  17                24  25  26    %

                        lrud = np.array([grid[UD-1,NS+1,EW], grid[UD,NS+1,EW], grid[UD+1,NS+1,EW],
                                         grid[UD-1,NS  ,EW],                   grid[UD+1,NS  ,EW],
                                         grid[UD-1,NS-1,EW], grid[UD,NS-1,EW], grid[UD+1,NS-1,EW],
                                         grid[UD-1,NS+1,EW+1], grid[UD,NS+1,EW+1], grid[UD+1,NS+1,EW+1],
                                         grid[UD-1,NS  ,EW+1], grid[UD,NS  ,EW+1], grid[UD+1,NS  ,EW+1],
                                         grid[UD-1,NS-1,EW+1], grid[UD,NS-1,EW+1], grid[UD+1,NS-1,EW+1],
                                         grid[UD-1,NS+1,EW-1], grid[UD,NS+1,EW-1], grid[UD+1,NS+1,EW-1],
                                         grid[UD-1,NS  ,EW-1], grid[UD,NS  ,EW-1], grid[UD+1,NS  ,EW-1],
                                         grid[UD-1,NS-1,EW-1], grid[UD,NS-1,EW-1], grid[UD+1,NS-1,EW-1]])

                        # count how many cell have data (i.e. values gt 0) & what their indices are
                        # -------------------------------------------------------------------------
                        ctaW = len(lrud[lrud > 0])
                        w_num = np.where(lrud > 0)

                        # if adjoining cells are unassigned, give cell value = objnum and increase objnum by one
                        if ctaW == 0:
                            grid[UD,NS,EW]=objnum
                            objnum=objnum+1
                        # if adjoining cells have values give the cell the min objnum
                        else:
                            grid[UD,NS,EW]=min(lrud[w_num])

        # START HERE ------------------------------------------------------------------------------------
        # It could have happened, that one storm has different objnums but we want
        # all cells in one storm to have the same value (numbername), so we go through the grid
        # again and look if adjoining non empty cells have different values.

        # this is evaluating from bottom to top (UpDown)
        for UD in range(0,tam2[0]-1):
            
            # this is evaluating from North to South  (lats-rows)
            for NS in range(0,tam2[1]-1):
                
                # this is evaluating from East to West (lons-cols)
                for EW in range(0,tam2[2]-1):

                    # print('UD=',UD,' NS=',NS,' EW=',EW,' grid[UD,NS,EW]=',grid[UD,NS,EW])
         
                    if grid[UD,NS,EW] != fillValue:
                        lrud2 = np.array([grid[UD-1,NS+1,EW], grid[UD,NS+1,EW], grid[UD+1,NS+1,EW],
                                          grid[UD-1,NS  ,EW], grid[UD,NS,EW] ,  grid[UD+1,NS  ,EW],
                                          grid[UD-1,NS-1,EW], grid[UD,NS-1,EW], grid[UD+1,NS-1,EW],
                                          grid[UD-1,NS+1,EW+1], grid[UD,NS+1,EW+1], grid[UD+1,NS+1,EW+1],
                                          grid[UD-1,NS  ,EW+1], grid[UD,NS  ,EW+1], grid[UD+1,NS  ,EW+1],
                                          grid[UD-1,NS-1,EW+1], grid[UD,NS-1,EW+1], grid[UD+1,NS-1,EW+1],
                                          grid[UD-1,NS+1,EW-1], grid[UD,NS+1,EW-1], grid[UD+1,NS+1,EW-1],
                                          grid[UD-1,NS  ,EW-1], grid[UD,NS  ,EW-1], grid[UD+1,NS  ,EW-1],
                                          grid[UD-1,NS-1,EW-1], grid[UD,NS-1,EW-1], grid[UD+1,NS-1,EW-1]])

                        # count how many cell have data (i.e. values gt 0) & what their indices are
                        # -------------------------------------------------------------------------
                        ctaW = len(lrud2[lrud2 > 0])
                        w_num = np.where(lrud2 > 0)

                        if ctaW != 0: 
                            # find if cells within the same storm have different objnum (std~=0)
                            if np.std(lrud2[w_num]) != 0.:
                                a = min(lrud2[w_num])
                                b = max(lrud2[w_num])
                                # set all cells with the different value to the minimum
                                grid[grid == b] = a
                        else:
                            print('findStormNew: There are no storms in this array')

        # Remove padding (1 row,col,lev to either side of grid) 
        grid_storm = np.array([grid[1:tam2[0]-1,1:tam2[1]-1,1:tam2[2]-1]])

        # make arrays with the IDs of all unique storms (id_storm) and counts for each one (npix_str)
        id_storm, npix_str = np.unique(grid_storm[grid_storm >= 0],return_counts=True)
        
    else:
        print('findStormNew: There are no reflectivities > 0 in this array')

    return id_storm, npix_str, grid_storm
