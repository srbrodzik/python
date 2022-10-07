pro findStorm,refl_3D=refl_3D,id_storm=id_storm,npix_str=npix_str,grid_storm=grid_storm

;*********************************************************************************
;here I am going to find in the 3D reflectivity matrix contiguous pixels with data
;matching the thereshold for reflectivty . Identify individual storms!!

;here I select only pixels with values greater than 0
tam=size(refl_3D)
cswath=lonarr(tam[1],tam[2],tam[3])
cswath[*,*,*]=-999l
donde=where(refl_3D gt 0.,cta)
if cta ne 0 then cswath[donde]=0l ;else stop ;here stops if there is no data!!!

;here I am going to expand the refl3D matrix 1 row,col,lev to either side
grid=lonarr(tam[1]+2,tam[2]+2,tam[3]+2)
grid[*,*,*]=-999l               ;value matching missing data
grid[1:tam[1],1:tam[2],1:tam[3]]=cswath

tam2=size(grid)
objnum=1l                        ;changes when a new storm is found

for UD=0,tam2[3]-1 do begin          ;;this is evaluating from bottom to top (UpDown)
   for NS=0,tam2[2]-1 do begin       ;; this is evaluating from North to South  (lats-rows)
      for EW=0,tam2[1]-1 do begin    ;; this is evaluating from East to West (lons-cols)
         
         if grid[EW,NS,UD] ne -999l then begin ;print,grid[EW,NS,UD]
            ;grab the cells around the cell of interest in all directions
            ;;NUMBERING OF SURROUNDING CELLS %%%
            ;    samelevel =  1   2   3    uplevel =  9  10  11    downlevel = 18  19  20    %
            ;                 4       5              12  13  14                21  22  23    %
            ;                 6   7   8              15  16  17                24  25  26    %
            lrud=[grid[EW-1,NS+1,UD], grid[EW,NS+1,UD], grid[EW+1,NS+1,UD], $
                  grid[EW-1,NS  ,UD],                   grid[EW+1,NS  ,UD], $
                  grid[EW-1,NS-1,UD], grid[EW,NS-1,UD], grid[EW+1,NS-1,UD], $
                  grid[EW-1,NS+1,UD+1], grid[EW,NS+1,UD+1], grid[EW+1,NS+1,UD+1], $
                  grid[EW-1,NS  ,UD+1], grid[EW,NS  ,UD+1], grid[EW+1,NS  ,UD+1], $
                  grid[EW-1,NS-1,UD+1], grid[EW,NS-1,UD+1], grid[EW+1,NS-1,UD+1], $
                  grid[EW-1,NS+1,UD-1], grid[EW,NS+1,UD-1], grid[EW+1,NS+1,UD-1], $
                  grid[EW-1,NS  ,UD-1], grid[EW,NS  ,UD-1], grid[EW+1,NS  ,UD-1], $
                  grid[EW-1,NS-1,UD-1], grid[EW,NS-1,UD-1], grid[EW+1,NS-1,UD-1]]

            w_num=where(lrud gt 0l,ctaW) ;count how many cell have data
 
            if ctaW eq 0l then begin ;if adjoining cells are unassigned, give cell value = objnum
               grid[EW,NS,UD]=objnum
               objnum=objnum+1l  ;increase the number to count for the next contiguous
            endif else begin    ;if adjoining cells have values give the cell the min objnum
               grid[EW,NS,UD]=min(lrud[w_num])
            endelse             ;if length...
         endif                  ; if(isfinite...

      endfor
   endfor
endfor

;%It could have happened, that one storms has different objnums but we want
;%all cells in one storm to have the same value (numbername), so we go through the grid
;%again and look if adjoining non empty cells have different values.

for UD=0,tam2[3]-1 do begin          ;;this is evaluating from bottom to top (UpDown)
   for NS=0,tam2[2]-1 do begin       ;; this is evaluating from North to South  (lats-rows)
      for EW=0,tam2[1]-1 do begin    ;; this is evaluating from East to West (lons-cols)

         if grid[EW,NS,UD] ne -999l then begin ;print,grid[EW,NS,UD]
            lrud2=[grid[EW-1,NS+1,UD], grid[EW,NS+1,UD], grid[EW+1,NS+1,UD], $
                   grid[EW-1,NS  ,UD],  grid[EW,NS,UD] ,  grid[EW+1,NS  ,UD], $
                   grid[EW-1,NS-1,UD], grid[EW,NS-1,UD], grid[EW+1,NS-1,UD], $
                   grid[EW-1,NS+1,UD+1], grid[EW,NS+1,UD+1], grid[EW+1,NS+1,UD+1], $
                   grid[EW-1,NS  ,UD+1], grid[EW,NS  ,UD+1], grid[EW+1,NS  ,UD+1], $
                   grid[EW-1,NS-1,UD+1], grid[EW,NS-1,UD+1], grid[EW+1,NS-1,UD+1], $
                   grid[EW-1,NS+1,UD-1], grid[EW,NS+1,UD-1], grid[EW+1,NS+1,UD-1], $
                   grid[EW-1,NS  ,UD-1], grid[EW,NS  ,UD-1], grid[EW+1,NS  ,UD-1], $
                   grid[EW-1,NS-1,UD-1], grid[EW,NS-1,UD-1], grid[EW+1,NS-1,UD-1]]

            w_num=where(lrud2 gt 0l,ctaW)
            if ctaW eq 0l then stop ;here is a problem!!!
            ;now try to find if cells within the same storm have different objnum (std~=0)
            if stddev(float(lrud2[w_num])) ne 0. then begin
               a=min(lrud2[w_num])
               b=max(lrud2[w_num])
               grid[where(grid eq b)]=a ;set all cells with the different value to the minimum
            endif
         endif                  ; if(isfinite...
      endfor
   endfor
endfor
grid_storm=(grid[1:tam2[1]-2,1:tam2[2]-2,1:tam2[3]-2])

;make a list with the IDs of all individual storms
ids_s=grid_storm[where(grid_storm ge 0l)]
sortedIDs=ids_s[sort(ids_s)]
pos=uniq(sortedIDs)
id_storm=sortedIDs[pos]

npix_str=lonarr(n_elements(id_storm))

for n=0,n_elements(id_storm)-1 do begin & $
  donde=where(grid_storm eq id_storm[n],cta) & $ ;;%count the total number of pixels in each storm
  npix_str[n]=cta  & $
endfor
;*********************************************************************************

end
