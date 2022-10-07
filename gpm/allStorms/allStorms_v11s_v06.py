#!/usr/bin/python3

import numpy as np
from findStormNew import findStorm
from lisotd2_v2 import listotd2_v2
from utilities import deg2km
from utilities import read_topo

def allStorms_v11s_v06(region,limits,path_in,path_out,topo_file,type):

    '''
    HISTORY
    Original code in 
       /home/disk/shear2/brodzik/matlab/ULLI.stormStats
    This code based on newest version of that IDL code:
       /home/disk/shear2/brodzik/IDL/gpm/allStorms/LandRuns/run_allStorms_wmp_v11s_v06.pro
    
    First uses findStorm code to locate all contiguous pixels with refl > 0
    Then uses findStorm when applying theresholds to account for the 
    classification (e.g., convective pixels and reflectivity >= 40 dBZ)
  
    NOTE from M Zuluaga about volume calculations (e.g. vol_Rain_All with units of 10^6 kg/s):
    "That unit of rainfall is called Rainfall Productivity that is [i.e., rain rate (kg m^-2 s^-1) x area (m2 )].
    This is to follow Romatschke and Houze 2011. If you work out the units, you'll get that.
    Multiply the rainfall rate (mm/hr) times the density of water (to get a volume)."
    '''

    # Define vars
    nlevels = 176
    delta_z = 0.125   # in km
  
    thr_aST=50000.      # threshold area to define Broad Stratiform regions
    thr_dbZ=40.         # reflectivity threshold for id the convective subset
    thr_hCV=10.         # height of the echo top for id deep cores
    thr_aCV=1000.       # area  of the echo core for id wide cores
    thr_dCV=5.          # min depth (maxHt-minHt) to be considered for DCC category
    
    STRA = 10
    CONV = 20
    OTHER = 30
    SHIS = 40

    pixKm = 5.5        ;;km
    pixArea = pixKm*pixKm
    npix_aST = fix(thr_aST/pixArea)
    pixDeg = 0.05      ;; deg
    secsPerHr = 3600.  ;; seconds

    # set this flag to 0/1 if you DO NOT/DO want masks added to input files
    makeCoreMasks = 1
  
    # read the topography data
    (topo_lat,topo_lon,DEM) = read_topo(topo_file)

    #variable definition for computation of Rain rates using Z-R method. Values from Romatscke and Houze 2010 
    aRR=140   # %Z-R relation values for a and b (Z=a*R^b)  other
    bRR=1.6   #
    aRRc=100  # %Z-R relation values for a and b (Z=a*R^b) convective
    bRRc=1.7  
    aRRs=200  # %Z-R relation values for a and b (Z=a*R^b) stratiform
    bRRs=1.49 

    #*************************NO NEED TO EDIT ANYTHING BELOW THIS LINE*************************
  
    #*******************************************************************************************
    #this are definition vector to store shape, rain parameters for each core and full storm
    info_DC=[]
    info_WC=[]
    info_DW=[]
    info_BS=[]
    info_SH=[]

    shape_Core_DC=fltarr(9,1) & shape_Full_DC=fltarr(9,1)
    shape_Core_WC=fltarr(9,1) & shape_Full_WC=fltarr(9,1)
    shape_Core_DW=fltarr(9,1) & shape_Full_DW=fltarr(9,1)
    shape_Core_BS=fltarr(9,1) & shape_Full_BS=fltarr(9,1)
    shape_Core_SH=fltarr(9,1) & shape_Full_SH=fltarr(9,1)

    rain_Core_DC=fltarr(7,1)  & rain_Full_DC=fltarr(7,1)
    rain_Core_WC=fltarr(7,1)  & rain_Full_WC=fltarr(7,1)
    rain_Core_DW=fltarr(7,1)  & rain_Full_DW=fltarr(7,1)
    rain_Core_BS=fltarr(7,1)  & rain_Full_BS=fltarr(7,1)
    rain_Core_SH=fltarr(7,1)  & rain_Full_SH=fltarr(7,1)

    rainTypeCore_DC=fltarr(6,1)  & rainTypeFull_DC=fltarr(6,1)
    rainTypeCore_WC=fltarr(6,1)  & rainTypeFull_WC=fltarr(6,1)
    rainTypeCore_DW=fltarr(6,1)  & rainTypeFull_DW=fltarr(6,1)
    rainTypeCore_BS=fltarr(6,1)  & rainTypeFull_BS=fltarr(6,1)
    rainTypeCore_SH=fltarr(6,1)  & rainTypeFull_SH=fltarr(6,1)

  ;;*******************************************************************************************
  ;;this are definition of matrices to count number of events in coarse grid
  res=0.5
  lonsC=findgen(fix(1l+(limits[3]-limits[1])/res))*res+limits[1]
  latsC=findgen(fix(1l+(limits[2]-limits[0])/res))*res+limits[0]
  nlonsC=fix(1l+(limits[3]-limits[1])/res)
  nlatsC=fix(1l+(limits[2]-limits[0])/res)

  freq_Full=lonarr(nlonsC,nlatsC,5)  ;;5 type of systems (DC,WC,DW,BS,SH)
  freq_Core=lonarr(nlonsC,nlatsC,5)  ;;5 type of systems (DC,WC,DW,BS,SH)

  ;;**************************************************************************
  ;;this is definition of matrices to calculate rain rates in fine grid
  res_f=0.05d
  lonsF=findgen(fix(1l+(limits[3]-limits[1])/res_f))*res_f+limits[1]
  latsF=findgen(fix(1l+(limits[2]-limits[0])/res_f))*res_f+limits[0]
  nlonsF=fix(1l+(limits[3]-limits[1])/res_f)
  nlatsF=fix(1l+(limits[2]-limits[0])/res_f)

  ;;*******************************************************************************************
  ;;Cummulative rainfal rate ;;***********************************************
  ;;rain_R11Full=fltarr(nlonsF,nlatsF,5)    ;;5 type of systems (DC,WC,DW,BS,SH)   ;; SRB
  ;;nRai_R11Full=intarr(nlonsF,nlatsF,5)    ;;5 type of systems (DC,WC,DW,BS,SH)   ;; SRB
  ;;rain_R11Core=fltarr(nlonsF,nlatsF,5)                                           ;; SRB
  ;;nRai_R11Core=intarr(nlonsF,nlatsF,5)                                           ;; SRB

  rain_NSRFull=fltarr(nlonsF,nlatsF,5)    ;;5 type of systems (DC,WC,DW,BS,SH)
  nRai_NSRFull=intarr(nlonsF,nlatsF,5)    ;;5 type of systems (DC,WC,DW,BS,SH)
  rain_NSRCore=fltarr(nlonsF,nlatsF,5)      
  nRai_NSRCore=intarr(nlonsF,nlatsF,5)      

  ;;rainCore_DC_R11=fltarr(6,1)   &   rainFull_DC_R11=fltarr(6,1)   ;; SRB
  ;;rainCore_WC_R11=fltarr(6,1)   &   rainFull_WC_R11=fltarr(6,1)   ;; SRB
  ;;rainCore_DW_R11=fltarr(6,1)   &   rainFull_DW_R11=fltarr(6,1)   ;; SRB
  ;;rainCore_BS_R11=fltarr(6,1)   &   rainFull_BS_R11=fltarr(6,1)   ;; SRB
  ;;rainCore_SH_R11=fltarr(6,1)   &   rainFull_SH_R11=fltarr(6,1)   ;; SRB

  rainCore_DC_NSR=fltarr(6,1)   &   rainFull_DC_NSR=fltarr(6,1)
  rainCore_WC_NSR=fltarr(6,1)   &   rainFull_WC_NSR=fltarr(6,1)
  rainCore_DW_NSR=fltarr(6,1)   &   rainFull_DW_NSR=fltarr(6,1)
  rainCore_BS_NSR=fltarr(6,1)   &   rainFull_BS_NSR=fltarr(6,1)
  rainCore_SH_NSR=fltarr(6,1)   &   rainFull_SH_NSR=fltarr(6,1)

  ;;*******************************************************************************************
  ;;this is definition of CFAD matrix
  n_refls=81
  alts_CFAD=findgen(nlevels)*delta_z  ;;altitudes in km
  refl_CFAD=findgen(n_refls)
  
  CFAD_Full=lonarr(n_refls,nlevels,5)  ;;5 type of systems (DC,WC,DW,BS,SH)
  CFAD_Core=lonarr(n_refls,nlevels,5)  ;;5 type of systems (DC,WC,DW,BS,SH)

  tmp=str_sep(type,'_')
  years=tmp[1] & meses=tmp[0] & region=tmp[2]

  path=path_in+years+'/'+meses
  cd,path
  files=findfile('*.nc',count=countNfiles)

  for ff=0l,countNfiles-1 do begin         ;; problem is at ff=1077, ss=84 and ssCV=5
     csa=str_sep(files[ff],'_')
     suffix = str_sep(csa[6],'.')
     filen=csa[0]+'_'+csa[1]+'_'+csa[2]+'_'+csa[3]+'_'+csa[4]+'_'+csa[5]+'_'+suffix[0]
     orbit=csa[5] & datetime=csa[2]
     print,'analyzing '+region+' region with allStorms_v11s_v06 for orbit='+orbit+' datetime='+datetime+$
           ' for file='+strtrim(ff+1,2)+'/'+strtrim(countNfiles,2)

     ;; Don't think we need this any more
     ;;n_col=0l & n_row=0l           ;;this comes from matlab, so col-fil are switched

     ;; MASK MOD requires file opened with write permission if makeCoreMasks is set
     if makeCoreMasks then begin
        ncid = ncdf_open(path+'/'+filen+'.nc',/write)
     endif else begin
        ncid = ncdf_open(path+'/'+filen+'.nc',/nowrite)
     endelse
        
     ;; get file dimensions
     lonDimID = ncdf_dimid(ncid,'lon')
     ncdf_diminq,ncid,lonDimID,name,nlons  ;; nlon is old n_col
     latDimID = ncdf_dimid(ncid,'lat')
     ncdf_diminq,ncid,latDimID,name,nlats  ;; nlat is old n_row
     altDimID = ncdf_dimid(ncid,'alt')
     ncdf_diminq,ncid,altDimID,name,nalts
     timDimID = ncdf_dimid(ncid,'time')
     ncdf_diminq,ncid,timDimID,name,ntimes

     ;; get vars
     latID = ncdf_varid(ncid,'lat')
     ncdf_varget,ncid,latID,lats
     lonID = ncdf_varid(ncid,'lon')
     ncdf_varget,ncid,lonID,lons
     hgts=findgen(nlevels)*delta_z           ;;this is kilometers CHANGE hgts to alts

     raintypeID = ncdf_varid(ncid,'rain_type_raw')        ;; CHANGE var name used to be rain_type_orig
     ncdf_varget,ncid,raintypeID,raintype
     ncdf_attget,ncid,raintypeID,'_FillValue',raintype_fillValue
     ncdf_attget,ncid,raintypeID,'no_rain_value',raintype_noRainValue
     ;;raintype_uwID = ncdf_varid(ncid,'rain_type_uw')      ;; 
     ;;ncdf_varget,ncid,raintype_uwID,raintype_uw
     shallow_rtypeID = ncdf_varid(ncid,'shallow_rain_type')
     ncdf_varget,ncid,shallow_rtypeID,shallow_raintype
     ;;apply info from raintype_uw (change some CONV to STRA where appropriate)
     ;;d_strt_uw=where(raintype_uw eq 1, cta_strt_uw)
     ;;if cta_strt_uw ne 0 then raintype[d_strt_uw]=10000000.
     
     surf_rainID = ncdf_varid(ncid,'near_surf_rain')
     ncdf_varget,ncid,surf_rainID,SrfRain    ;; CHANGE var name to near_surf_rain
     ncdf_attget,ncid,surf_rainID,'_FillValue',SrfRain_fillValue
     
     reflID = ncdf_varid(ncid,'refl')        ;; CHANGE var name used to be maxdz
     ncdf_varget,ncid,reflID,refl_3D
     ncdf_attget,ncid,reflID,'_FillValue',refl_3D_fillValue

     ;;this uses the native raintype data to select more categories..
     d_strt=where(raintype ge 10000000 and raintype lt 20000000,cta_strt)
     d_conv=where(raintype ge 20000000 and raintype lt 30000000,cta_conv)
     d_NoShal=where(shallow_raintype eq 0,cta_NoShal)                             ;;no shallow rain (0)
     d_Shal=where(shallow_raintype eq 10 or shallow_raintype eq 11 or $
                  shallow_raintype eq 20 or shallow_raintype eq 21,cta_Shal)      ;;all shallow rain (10,11,20,21)
     d_ShIs=where(shallow_raintype eq 10 or shallow_raintype eq 11,cta_ShIs)      ;;all shallow Isolated rain (10,11)
     d_ShNi=where(shallow_raintype eq 20 or shallow_raintype eq 21,cta_ShNi)      ;;all shallow Non-Isolated rain (20,21)
     d_othe=where(raintype ge 30000000,cta_othe)
  
     d_noRa=where(raintype eq raintype_noRainValue,cta_noRa)  ;;generalized
     d_miss=where(raintype eq raintype_fillValue,cta_miss)    ;;generalized
     
     if cta_strt ne 0 then raintype[d_strt]=STRA  ;;stratiform
     if cta_conv ne 0 then raintype[d_conv]=CONV  ;;convective
     if cta_othe ne 0 then raintype[d_othe]=OTHER ;;Others

     ;;this are the definitions for shallow rain
     ;;if cta_Shal ne 0 then raintype[d_Shal]=SHIS   ;;all is Shallow
     if cta_ShIs ne 0 then raintype[d_ShIs]=SHIS     ;;all is Shallow Isolated
     ;;if cta_ShNi ne 0 then raintype[d_ShNi]=SHIS   ;;all is Shallow Non-Isolated

     rain_type3D=lonarr(nlons,nlats,nlevels)
     for i=0,nlevels-1 do rain_type3D[*,*,i]=raintype

     hgts_3D=fltarr(nlons,nlats,nlevels)
     for i=0l,nlevels-1l do hgts_3D[*,*,i]=hgts[i]

     ;; define new variables for masks if makeCoreMasks is set
     if makeCoreMasks then begin
        ncdf_control,ncid,/redef
        
        mask_missing_value = -99.
        result = ncdf_varid(ncid,'bsr_mask_str')
        if result eq -1 then begin
           bsr_id = ncdf_vardef(ncid,'bsr_mask_str',[lonDimID,latDimID,timDimID],/float,gzip=5)
           ncdf_attput,ncid,bsr_id,'_FillValue',mask_missing_value,/float
           ncdf_attput,ncid,bsr_id,'units','none',/char
           ncdf_attput,ncid,bsr_id,'long_name','BroadStratiform Core mask for strong thresholds',/char
        endif 
        result = ncdf_varid(ncid,'dcc_mask_str')
        if result eq -1 then begin
           dcc_id = ncdf_vardef(ncid,'dcc_mask_str',[lonDimID,latDimID,timDimID],/float,gzip=5)
           ncdf_attput,ncid,dcc_id,'_FillValue',mask_missing_value,/float
           ncdf_attput,ncid,dcc_id,'units','none',/char
           ncdf_attput,ncid,dcc_id,'long_name','DeepConvective Core mask for strong thresholds',/char
        endif 
        result = ncdf_varid(ncid,'dwc_mask_str')
        if result eq -1 then begin
           dwc_id = ncdf_vardef(ncid,'dwc_mask_str',[lonDimID,latDimID,timDimID],/float,gzip=5)
           ncdf_attput,ncid,dwc_id,'_FillValue',mask_missing_value,/float
           ncdf_attput,ncid,dwc_id,'units','none',/char
           ncdf_attput,ncid,dwc_id,'long_name','DeepWideConvective Core mask for strong thresholds',/char
        endif 
        result = ncdf_varid(ncid,'wcc_mask_str')
        if result eq -1 then begin
           wcc_id = ncdf_vardef(ncid,'wcc_mask_str',[lonDimID,latDimID,timDimID],/float,gzip=5)
           ncdf_attput,ncid,wcc_id,'_FillValue',mask_missing_value,/float
           ncdf_attput,ncid,wcc_id,'units','none',/char
           ncdf_attput,ncid,wcc_id,'long_name','WideConvective Core mask for strong thresholds',/char
        endif 
        ;;result = ncdf_varid(ncid,'shi_mask_str')
        ;;if result eq -1 then begin
        ;;   shi_id = ncdf_vardef(ncid,'shi_mask',[lonDimID,latDimID,timDimID],/float,gzip=5)
        ;;   ncdf_attput,ncid,shi_id,'_FillValue',mask_missing_value,/float
        ;;   ncdf_attput,ncid,shi_id,'units','none',/char
        ;;   ncdf_attput,ncid,shi_id,'long_name','ShallowIsolated Core mask',/char
        ;;endif
        result = ncdf_varid(ncid,'storm_mask_str')
        if result eq -1 then begin
           storm_id = ncdf_vardef(ncid,'storm_mask_str',[lonDimID,latDimID,timDimID],/float,gzip=5)
           ncdf_attput,ncid,storm_id,'_FillValue',mask_missing_value,/float
           ncdf_attput,ncid,storm_id,'units','none',/char
           ncdf_attput,ncid,storm_id,'long_name','Storm mask for strong thresholds',/char
        endif 
           
        ;; add global attribute descriptors for masks
        ncdf_attput,ncid,/global,'BroadStratiform_Criteria_Strong','contiguous stratiform echo >= 50,000km^2',/char
        ncdf_attput,ncid,/global,'DeepConvective_Criteria_Strong','contiguous, convective, 40dBZ echos with max height >= 10km',/char
        ncdf_attput,ncid,/global,'WideConvective_Criteria_Strong','contiguous, convective, 40dBZ echos with max horizontal dim >= 1000km^2',/char
        ncdf_attput,ncid,/global,'DeepWideConvective_Criteria_Strong','meets both DeepConvective_Criteria_Strong and WideConvective_Criteria_Strong',/char
        
        ;; put input file in data mode
        ncdf_control,ncid,/endef

        ;; allocate space for mask arrays
        bsr_mask = make_array(nlons,nlats,ntimes,/float,value=mask_missing_value)
        dcc_mask = make_array(nlons,nlats,ntimes,/float,value=mask_missing_value)
        dwc_mask = make_array(nlons,nlats,ntimes,/float,value=mask_missing_value)
        wcc_mask = make_array(nlons,nlats,ntimes,/float,value=mask_missing_value)
        ;;shi_mask = make_array(nlons,nlats,ntimes,/float,value=mask_missing_value)
        storm_mask = make_array(nlons,nlats,ntimes,/float,value=mask_missing_value)
     endif 
        
     ;; accumulate counts of extreme events
     num_bsr = 0
     num_dcc = 0
     num_dwc = 0
     num_wcc = 0
     ;;num_shi = 0
     num_storm = 0
     
     ;;Start Classification 
     ;;*********************************************************************************************
     refl_Rain=refl_3D
     where_raining=where(refl_Rain ge 0.,cta_0) ;;check for all pixels with >=0. dBZ
      
     if cta_0 gt 2 then begin             ;;here analyze only volumes with reflect >=0.
        ;;%%%%we call the function findStorm.m which searches for unique storms in the  swath
        refl_Rain[where(refl_3D lt 0.)]=refl_3D_fillValue ;;% threshold for reflectivitity ID everthing dBZ>0
        ;; -------------------- SRB Check why including missing val causes errors -----------------------
        findStormNew,refl_3D=refl_Rain,$
                     ;;refl_3D_fillValue=refl_3D_fillValue,$
                     id_storm=id_storm,$
                     npix_str=npix_str,$
                     grid_storm=grid_storm
        
        searchNaN=where(grid_storm lt 0, nanCnt)
        if nanCnt gt 0 then grid_storm[searchNaN]=0l    ;;%set NaNs to zero in the  grid matrix 

        ;;*******************
        ;;*******************
        print,'Id-Stratiform'
        ;;*******************
        ;;*******************
        ;;Identify The Broad Stratiform regions
        ;;identify only volumes with more than npix_aST (~1650 pixels) (~50000.km2/(5.5km*5.5km) Broad Stratiform class 
        donde_BrdStr=where(npix_str gt 1000l,ctaBrdStr) ;;1000 just to be conservative
        ;; SRB NEW VERSION OF THESE LINES
        ;;identify only volumes with more than npix_aST (~1650 pixels) (~50000.km2/(5.5km*5.5km) Broad Stratiform class 
        ;;donde_BrdStr=where(npix_str gt 1400l,ctaBrdStr) ;;1400 just to be conservative

        ;;here it is going to identify Stratiform subset
        if cta_strt ne 0 and ctaBrdStr gt 0 then begin ;;here analyze only stratiform and >=0.
 
           for ss=0l,ctaBrdStr-1 do begin  ;;only goes thru the volumes having more than 1000 pixels
              s_idF=id_storm[donde_BrdStr[ss]]
              w_idF=where(grid_storm eq s_idF,cta1)
              ;; make sure size of w_idF is same as npix_str for that storm id
              if cta1 ne npix_str[donde_BrdStr[ss]] then stop  ;; just to check!
              ;; create new array with storm id at appropriate indices
              singlestormgrid=lonarr(nlons,nlats,nlevels)
              singlestormgrid[w_idF]=s_idF
      
              ;;Here I will subset the singleFULL storm found to reduce the size of the matrix
              ;; find lon indices of storm in singlestormgrid
              total_lonFull=total(total(singlestormgrid,2),2) 
              d_lonsFull=where(total_lonFull gt 0l,nlonsFull) ;;this are the horizontal positions within a selcted contiguous area
              ;; find lat indices of storm in singlestormgrid
              total_latFull=total(total(singlestormgrid,1),2) 
              d_latsFull=where(total_latFull gt 0l,nlatsFull)

              ;; create array of actual storm dims containing singlestormgrid values
              singlestormgrid_Full=lonarr(nlonsFull,nlatsFull,nlevels)
              singlestormgrid_Full=singlestormgrid[d_lonsFull,d_latsFull,*]
              w_idF=where(singlestormgrid_Full eq s_idF,cta_Full)
              ;; these are actual lons and lats for singlestormgrid_Full array
              lonsFull_sub=lons[d_lonsFull]
              latsFull_sub=lats[d_latsFull]

              undefine,singlestormgrid

              ;;Identify the Stratiform pixels in storm (refl >= 0 and All stratiform pixels)
              singlestormgridStratiform=singlestormgrid_Full
              singlestormgridStratiform[where(refl_3D[d_lonsFull,d_latsFull,*] lt 0.)]=0l         ;;all dbZ greater than 0.
              singlestormgridStratiform[where(rain_type3D[d_lonsFull,d_latsFull,*] ne STRA)]=0l   ;;all stratiform pixels=STRA
               
              ;;2D horizontal array of number of Stratiform pixels within storm
              grid_sumST=total(singlestormgridStratiform,3) 
              dondeST=where(grid_sumST gt 0,pixelsumST) ;;pixelsum=number of pixels in a 2D proj

              if pixelsumST ge 2 then begin ;; look at events with 2 or more stratiform pixels
                 lonST=total(total(singlestormgridStratiform,2),2) 
                 lonC_ST=(max(lonsFull_sub[where(lonST gt 0l)])+min(lonsFull_sub[where(lonST gt 0l)]))/2. ;;center of Stratiform

                 latST=total(total(singlestormgridStratiform,1),2) 
                 latC_ST=(max(latsFull_sub[where(latST gt 0l)])+min(latsFull_sub[where(latST gt 0l)]))/2. ;;center of Stratiform

                 size_pixels=deg2km(pixDeg,lonC_ST,latC_ST)
                 area_ST=pixelsumST*(size_pixels[0]*size_pixels[1]) ;;stratiform area in km  
              endif else area_ST=0.

              ;;if there are pixels within individual cluster accomplishing the Strat criteria
              if area_ST ge thr_aST then begin ;;
                 ;;Does the stratiform area belongs to a single storm? (acomplish contiguous pixel condition)
                 ;;Again, uses the stratiform pixels to identify contiguous pixels
                 findStormNew,refl_3D=singlestormgridStratiform,$
                              ;;refl_3D_fillValue=refl_3D_fillValue,$
                              id_storm=id_ST,$
                              npix_str=npix_ST,$
                              grid_storm=grid_ST
                 searchNaN_ST=where(grid_ST lt 0, nanCnt)
                 if nanCnt gt 0 then grid_ST[searchNaN_ST]=0l ;;%set NaNs to zero in the  grid matrix
                 
                 ;;identify only volumes with more than npix_aST (~1650 pixels) (~50000./(5.5*5.5) Broad Stratiform class 
                 donde_BrdStr2=where(npix_ST gt 1000l,ctaBrdStr2) ;;1000 just to be conservative
                 ;; SRB NEW VERSION OF THIS LINE
                 ;;donde_BrdStr2=where(npix_ST gt 1400l,ctaBrdStr2) ;;1400 just to be conservative

                 areaIsFULL=0. & lonCIsFULL=0. & latCIsFull=0. ;;this is to locate only ONE full storm

                 for ssST=0l,ctaBrdStr2-1 do begin ;;only goes thru volumes greater than 1000 pixels (speed the process)
                    s_idST=id_ST[donde_BrdStr2[ssST]]
                    w_idST=where(grid_ST eq s_idST,cta_ST)
                    ;; make sure size of w_idST is same as npix_ST for that storm id
                    if cta_ST ne npix_ST[donde_BrdStr2[ssST]] then stop  ;; just to check!

                    ;; create new array with storm id at appropriate indices
                    singlestormgrid_ST=lonarr(nlonsFull,nlatsFull,nlevels)
                    singlestormgrid_ST[w_idST]=s_idST
                    ;;2D horizontal array of number of Stratiform pixels within storm
                    grid_sum_ST_ST=total(singlestormgrid_ST,3) 
                    dondeST_ST=where(grid_sum_ST_ST gt 0,pixelsumST_ST) ;;pixelsum=number of pixels in a 2D proj

                    lonST=total(total(singlestormgrid_ST,2),2) 
                    lonC_ST=(max(lonsFull_sub[where(lonST gt 0l)])+min(lonsFull_sub[where(lonST gt 0l)]))/2. ;;center of Stratiform

                    latST=total(total(singlestormgrid_ST,1),2) 
                    latC_ST=(max(latsFull_sub[where(latST gt 0l)])+min(latsFull_sub[where(latST gt 0l)]))/2. ;;center of Stratiform

                    size_pixels=deg2km(pixDeg,lonC_ST,latC_ST)
                    area_BS=pixelsumST_ST*(size_pixels[0]*size_pixels[1]) ;;Broad Stratiform area in km  
         
                    ;;now it tries to identify the real Broad Stratiform contiguous area
                    if area_BS ge thr_aST then begin ;;

                       ;; increment num_bsr
                       num_bsr = num_bsr + 1

                       ;; add bsr core number to bsr_mask if makeCoreMasks is set
                       if makeCoreMasks then begin
                          bsr_mask_sub = bsr_mask[d_lonsFull,d_latsFull,0]
                          bsr_mask_sub[dondeST_ST]=num_bsr
                          bsr_mask[d_lonsFull,d_latsFull,0] = bsr_mask_sub
                          undefine,bsr_mask_sub
                       endif 
                       
                       ;;it calculate the dimensions for Broad Stratiform pixels in the cluster
                       ;;******************************************************************************************
                       dim_lonBS=(max(lonsFull_sub[where(lonST gt 0l)])-min(lonsFull_sub[where(lonST gt 0l)]))+pixDeg
                       dim_latBS=(max(latsFull_sub[where(latST gt 0l)])-min(latsFull_sub[where(latST gt 0l)]))+pixDeg

                       hgt_sumBS=total(total(singlestormgrid_ST,2),1)
                       dim_hgtBS=max(hgts[where(hgt_sumBS gt 0l)])-min(hgts[where(hgt_sumBS gt 0l)])
                       dim_topBS=max(hgts[where(hgt_sumBS gt 0l)])   ;;###############
                       dim_botBS=min(hgts[where(hgt_sumBS gt 0l)])   ;;###############

                       ;;calculates the elevation of the terrain for the center of the storm  !21
                       id_top1=(where(topo_lon ge lonC_ST))[0] & id_top2=(where(topo_lat ge latC_ST))[0]
                       terr_hgtBS=DEM[id_top1,id_top2]                            ;;elevation in meters
                       if terr_hgtBS eq 0 then land_oceanBS=0 else land_oceanBS=1 ;;ocean=0 or land=1!

                       tmp_shapeBS=[lonC_ST,latC_ST,area_BS,dim_topBS,dim_botBS,dim_lonBS,dim_latBS,terr_hgtBS,land_oceanBS]

                       ;;******************************************************************************************
                       ;;Compute the Shape parameters and info of full storm (all pixels in the cluster)
                       ;;******************************************************************************************
                       grid_sum=total(singlestormgrid_Full,3)    ;;2D horiz array of number of all pixels within storm
                       donde=where(grid_sum gt 0,pixelsum)       ;;pixelsum=number of pixels in a 2D proj

                       ;; increment num_storm
                       num_storm = num_storm + 1

                       ;; add storm to storm_mask if makeCoreMasks is set
                       if makeCoreMasks then begin
                          storm_mask_sub = storm_mask[d_lonsFull,d_latsFull,0]
                          storm_mask_sub[donde]=num_storm
                          storm_mask[d_lonsFull,d_latsFull,0] = storm_mask_sub
                          undefine,storm_mask_sub
                       endif 

                       lon_sum=total(total(singlestormgrid_Full,2),2)   
                       cen_lon=(max(lonsFull_sub[where(lon_sum gt 0l)])+min(lonsFull_sub[where(lon_sum gt 0l)]))/2.     ;;longitude coordinate
                       dim_lon=(max(lonsFull_sub[where(lon_sum gt 0l)])-min(lonsFull_sub[where(lon_sum gt 0l)]))+pixDeg ;;longitudinal extents

                       lat_sum=total(total(singlestormgrid_Full,1),2)   
                       cen_lat=(max(latsFull_sub[where(lat_sum gt 0l)])+min(latsFull_sub[where(lat_sum gt 0l)]))/2.     ;;latitude coordinate
                       dim_lat=(max(latsFull_sub[where(lat_sum gt 0l)])-min(latsFull_sub[where(lat_sum gt 0l)]))+pixDeg ;;latitudinal extents

                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       area=pixelsum*size_pixels[0]*size_pixels[1] ;;in km horizontal area of selected storm
     
                       hgt_sum=total(total(singlestormgrid_Full,2),1)
                       dim_hgt=(max(hgts[where(hgt_sum gt 0l)])-min(hgts[where(hgt_sum gt 0l)])) ;;%vertical extent in km
                       dim_top=max(hgts[where(hgt_sum gt 0l)])                                   ;;###############
                       dim_bot=min(hgts[where(hgt_sum gt 0l)])                                   ;;###############

                       ;;%calculates the elevation of the terrain for the center of the storm  !21
                       id_top1=(where(topo_lon ge cen_lon))[0] & id_top2=(where(topo_lat ge cen_lat))[0]
                       terr_hgt=DEM[id_top1,id_top2]                        ;;elevation in meters
                       if terr_hgt eq 0 then land_ocean=0 else land_ocean=1 ;;mask  ocean=0 or land=1!

                       tmp_shape=[cen_lon,cen_lat,area,dim_top,dim_bot,dim_lon,dim_lat,terr_hgt,land_ocean]

                       ;;**********************************************************************************************
                       ;;Statistics for rain types and rain rates Calculated using Near Surface Rain product
                       ;;Statistics over Broad Stratiform Area!
                       SrfRain_FULL=SrfRain[d_lonsFull,d_latsFull,*]
                       raintypeFULL=raintype[d_lonsFull,d_latsFull,*]
                       refl_3D_FULL=refl_3D[d_lonsFull,d_latsFull,*]
                       hgts_3D_FULL=hgts_3D[d_lonsFull,d_latsFull,*]

                       stratconv=raintypeFULL[dondeST_ST]                     ;;type of rain in each 2D pixel that compose the storm
                       strats=where(stratconv eq STRA,ctaStr)                 ;;stratiform
                       convec=where(stratconv eq CONV,ctaCon)                 ;;convective
                       others=where(stratconv ge OTHER,ctaOth)                ;;other type
                       noRain=where(stratconv eq raintype_noRainValue,ctaNoR) ;;no rain
                       missin=where(stratconv eq raintype_fillValue,ctaMis)   ;;missing value

                       rainBS=SrfRain_FULL[dondeST_ST] ;;this is based on Near Surface Rain
                       rain_nomiss=where(rainBS ne SrfRain_fillValue,Rmiss)
                       
                       ;;here I calculate moments of simple rain within storm
                       if Rmiss ge 2 then rain_momentBS=[mean(rainBS[rain_nomiss]),stdev(rainBS[rain_nomiss]),$
                                                         max(rainBS[rain_nomiss]),min(rainBS[rain_nomiss]),float(Rmiss),$
                                                         float(ctaStr),float(ctaCon)] $
                       else rain_momentBS=[-9999.,-9999.,-9999.,-9999.,-9999.,-9999.,-9999.]

                       ;;here I calculate rainrate sums in [mm/hr]
                       if ctaStr ne 0 then RTo_stra=total(rainBS[strats]) else RTo_stra=0.  ;;total stratiform rain
                       if ctaCon ne 0 then RTo_conv=total(rainBS[convec]) else RTo_conv=0.  ;;total convective rain
                       if ctaOth ne 0 then RTo_othe=total(rainBS[others]) else RTo_othe=0.  ;;total other rain
                       if ctaNoR ne 0 then RTo_noRa=total(rainBS[noRain]) else RTo_noRa=0.  ;;total no_rain rain
                       total_RainAll=RTo_stra+RTo_conv+RTo_othe+RTo_noRa
                       total_RainCSs=RTo_stra+RTo_conv

                       m_RainAll=total_RainAll/(ctaStr+ctaCon+ctaOth+ctaNoR+ctaMis)           ;;mean rainfall all pixels
                       if ctaStr ne 0 then m_RainStrt=RTo_stra/ctaStr else m_RainStrt=0.      ;;mean stratiform rain
                       if ctaCon ne 0 then m_RainConv=RTo_conv/ctaCon else m_RainConv=0.      ;;mean convective rain
                       if ctaStr ne 0 or ctaCon ne 0 then m_RainCSs=total_RainCSs/(ctaStr+ctaCon) $
                       else m_RainCSs=-999.

                       ;;here calculate volumetric values of rain  ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,lonC_ST,latC_ST)
                       vol_Rain_All=total_RainAll*(size_pixels[0]*size_pixels[1])/secsPerHr  
                       vol_Rain_Str=RTo_stra*(size_pixels[0]*size_pixels[1])/secsPerHr 
                       vol_Rain_Con=RTo_conv*(size_pixels[0]*size_pixels[1])/secsPerHr
   
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       statsRain_BS=[m_RainAll,m_RainStrt,m_RainConv,vol_Rain_All,vol_Rain_Str,vol_Rain_Con]
 
                       ;;******************************************************************************************
                       ;;statisitcs over the full region of the storm
                       stratconv=raintypeFULL[donde]                          ;;type of rain in each 2D pixel that compose the storm
                       strats=where(stratconv eq STRA,ctaStr)                 ;;stratiform
                       convec=where(stratconv eq CONV,ctaCon)                 ;;convective
                       others=where(stratconv ge OTHER,ctaOth)                ;;other type
                       noRain=where(stratconv eq raintype_noRainValue,ctaNoR) ;;no rain
                       missin=where(stratconv eq raintype_fillValue,ctaMis)   ;;missing value

                       ;;here I calculate moments of simple rain within storm
                       rain=SrfRain_FULL[donde] ;;this is based on Near Surface Rain
                       rain_nomiss=where(rain ne SrfRain_fillValue,Rmiss)
                       if Rmiss ge 2 then rain_moment=[mean(rain[rain_nomiss]),stdev(rain[rain_nomiss]),$
                                                       max(rain[rain_nomiss]),min(rain[rain_nomiss]),float(Rmiss),$
                                                       float(ctaStr),float(ctaCon)] $
                       else rain_moment=[-9999.,-9999.,-9999.,-9999.,-9999.,-9999.,-9999.]

                       ;;here I calculate rainrate sums in [mm/hr]
                       if ctaStr ne 0 then RTo_stra=total(rain[strats]) else RTo_stra=0.  ;;total stratiform rain
                       if ctaCon ne 0 then RTo_conv=total(rain[convec]) else RTo_conv=0.  ;;total convective rain
                       if ctaOth ne 0 then RTo_othe=total(rain[others]) else RTo_othe=0.  ;;total other rain
                       if ctaNoR ne 0 then RTo_noRa=total(rain[noRain]) else RTo_noRa=0.  ;;total no_rain rain
                       total_RainAll=RTo_stra+RTo_conv+RTo_othe+RTo_noRa
                       total_RainCSs=RTo_stra+RTo_conv

                       m_RainAll=total_RainAll/(ctaStr+ctaCon+ctaOth+ctaNoR+ctaMis)       ;;mean rainfall all pixels
                       if ctaStr ne 0 then m_RainStrt=RTo_stra/ctaStr else m_RainStrt=0.  ;;mean stratiform rain
                       if ctaCon ne 0 then m_RainConv=RTo_conv/ctaCon else m_RainConv=0.  ;;mean convective rain
                       if ctaStr ne 0 or ctaCon ne 0 then m_RainCSs=total_RainCSs/(ctaStr+ctaCon) else m_RainCSs=-999.

                       ;;here calculate volumetric values of rain  ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       vol_Rain_All=total_RainAll*(size_pixels[0]*size_pixels[1])/secsPerHr  
                       vol_Rain_Str=RTo_stra*(size_pixels[0]*size_pixels[1])/secsPerHr 
                       vol_Rain_Con=RTo_conv*(size_pixels[0]*size_pixels[1])/secsPerHr
   
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       statsRain=[m_RainAll,m_RainStrt,m_RainConv,vol_Rain_All,vol_Rain_Str,vol_Rain_Con]
            
                       ;;**********************************************************************************************
                       ;;Computation of rainfall accumulation per Full and Core storms
                       ;;**********************************************************************************************
                       grid_storm_FULL=lonarr(nlonsFull,nlatsFull,nlevels)
                       grid_storm_FULL=grid_storm[d_lonsFull,d_latsFull,*]

                       ;;----------------------------------------------------------------------------------------------
                       ;;MOVED THIS CODE INTO IF LOOP WHERE IT NEEDS TO BE TO GET RID OF NaN OUTPUTS
                       ;;Statistics for rain rates Calculated using Z-R method (Romatschke and Houze 2010)
                       ;;we need to get all the needed information for the lowest pixels in each storm
                       ;;locate the indices of the lowest pixels with data   (Full Storm) Romatsche and Houze 2011
                       ;;;;RTo_conv_R11=0. & ctaConv_R11=0l  ;;total convective rain  ;; SRB
                       ;;;;RTo_stra_R11=0. & ctaStra_R11=0l  ;;total stratiform rain  ;; SRB
                       ;;;;RTo_othe_R11=0. & ctaOthe_R11=0l  ;;total other rain       ;; SRB
                       ;;;;RTo_noRa_R11=0. & ctaNoRa_R11=0l  ;;total no_rain rain     ;; SRB

                       ;;Near Surface Rain
                       ;;RTo_conv_NSR=0. & ctaConv_NSR=0l ;;total convective rain
                       ;;RTo_stra_NSR=0. & ctaStra_NSR=0l ;;total stratiform rain
                       ;;RTo_othe_NSR=0. & ctaOthe_NSR=0l ;;total other rain
                       ;;RTo_noRa_NSR=0. & ctaNoRa_NSR=0l ;;total no_rain rain
                       ;;----------------------------------------------------------------------------------------------

                       ;;print,'BEFORE STORM CALCS: ss=',string(strtrim(ss,2)),' and ssST=',string(strtrim(ssST,2)), $
                       ;;      ' and area=',string(strtrim(area,2)),' and areaIsFULL=',string(strtrim(areaIsFULL,2))

                       ;;here also is accounted n_pixels in matrix matching with storm location and locate indices (lowest pixels
                       ;;this is only evaluated ONCE! (accouting for a single FULL storm)
                       if area ne areaIsFULL and cen_lon ne lonCIsFULL and cen_lat ne latCIsFull then begin

                          ;;print,'     RUNNING THROUGH CALCS FOR STORM - SHOULD BE ONLY ONCE'
                          
                          ;;Statistics for rain rates Calculated using Z-R method 
                          ;;RTo_conv_R11=0. & ctaConv_R11=0l  ;;total convective rain  ;; SRB
                          ;;RTo_stra_R11=0. & ctaStra_R11=0l  ;;total stratiform rain  ;; SRB
                          ;;RTo_othe_R11=0. & ctaOthe_R11=0l  ;;total other rain       ;; SRB
                          ;;RTo_noRa_R11=0. & ctaNoRa_R11=0l  ;;total no_rain rain     ;; SRB

                          ;;Statistics for rain rates Calculated using Near Surface Rain
                          RTo_conv_NSR=0. & ctaConv_NSR=0l ;;total convective rain
                          RTo_stra_NSR=0. & ctaStra_NSR=0l ;;total stratiform rain
                          RTo_othe_NSR=0. & ctaOthe_NSR=0l ;;total other rain
                          RTo_noRa_NSR=0. & ctaNoRa_NSR=0l ;;total no_rain rain
                          
                          for i=0l,pixelsum-1l do begin   
                             col=donde[i] mod nlonsFull                 ;;column ID of the pixel
                             fil=long(fix(donde[i]/float(nlonsFull)))   ;;row    ID of the pixel
                 
                             nearSrfR_Org=SrfRain_FULL[col,fil]

                             t_col=where(lonsC le lonsFull_sub[col],ctaC) ;;count the pixels
                             t_row=where(latsC le latsFull_sub[fil],ctaR) 
                             if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                colBS=(reverse(t_col))[0] 
                                rowBS=(reverse(t_row))[0] 
                                freq_Full[colBS,rowBS,3]=freq_Full[colBS,rowBS,3]+1l 

                                ;;locate indices (x,y,z location of individual pixels within the FULL storm)
                                ;;pila=reform(singlestormgrid_Full[col,fil,*])
                                pila=reform(grid_storm_FULL[col,fil,*])
                                w_CV=where(pila ne 0,ctaHgt)
                                if ctaHgt ge 1 then begin 
                                   ;;distance between lowest pixel and ground
                                   id_top1=(where(topo_lon ge lonsFull_sub[col]))[0]
                                   id_top2=(where(topo_lat ge latsFull_sub[fil]))[0]

                                   ;;%we sort out all the pixels that are higher than 2.5km above the ground
                                   ;;if (hgts[w_CV[0]]-float(DEM[id_top1,id_top2])/1000.) le 2.5 then begin    ;; SRB
                                      ;;locate the pixel in the nearest fine grid cell
                                      tmp_col=(where(float(lonsF) eq lonsFull_sub[col],ctaC))[0]
                                      tmp_row=(where(float(latsF) eq latsFull_sub[fil],ctaR))[0]

                                      if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                      ;;if ctaC eq 0 or ctaR eq 0 then stop
                          
                                         ;;if nearSrfR_Org ne -9999.00 and refl_3D_FULL[col,fil,w_CV[0]] ne -999.0 then begin               ;; SRB
                                         if nearSrfR_Org ne SrfRain_fillValue then begin
                                            ;;reflectivZ=10^(refl_3D_FULL[col,fil,w_CV[0]]*0.1)  ;;%convert from dBZ to Z                   ;; SRB

                                            rain_NSRFull[tmp_col,tmp_row,3]=rain_NSRFull[tmp_col,tmp_row,3]+nearSrfR_Org
                                            nRai_NSRFull[tmp_col,tmp_row,3]=nRai_NSRFull[tmp_col,tmp_row,3]+1l

                                            ;;here I create accumulated rain rate vectors
                                            if raintypeFULL[col,fil] eq CONV then begin ;;convective rain
                                               RTo_conv_NSR=RTo_conv_NSR+nearSrfR_Org                     &  ctaConv_NSR=ctaConv_NSR+1l
                                               ;;RTo_conv_R11=RTo_conv_R11+(reflectivZ/aRRc)^(1/bRRc)       &  ctaConv_R11=ctaConv_R11+1l   ;; SRB
                                            
                                               ;;rain_R11Full[tmp_col,tmp_row,3]=rain_R11Full[tmp_col,tmp_row,3]+(reflectivZ/aRRc)^(1/bRRc) ;; SRB
                                               ;;nRai_R11Full[tmp_col,tmp_row,3]=nRai_R11Full[tmp_col,tmp_row,3]+1l                         ;; SRB
                                            endif else begin
                                               if raintypeFULL[col,fil] eq STRA then begin ;;stratiform rain
                                                  RTo_stra_NSR=RTo_stra_NSR+nearSrfR_Org                  &  ctaStra_NSR=ctaStra_NSR+1l
                                                  ;;RTo_stra_R11=RTo_stra_R11+(reflectivZ/aRRs)^(1/bRRs)    &  ctaStra_R11=ctaStra_R11+1l   ;; SRB
                                               
                                                  ;;rain_R11Full[tmp_col,tmp_row,3]=rain_R11Full[tmp_col,tmp_row,3]+(reflectivZ/aRRs)^(1/bRRs)   ;; SRB 
                                                  ;;nRai_R11Full[tmp_col,tmp_row,3]=nRai_R11Full[tmp_col,tmp_row,3]+1l                      ;; SRB
                                               endif else begin
                                                  if raintypeFULL[col,fil] ge OTHER then begin ;;other rain
                                                     RTo_othe_NSR=RTo_othe_NSR+nearSrfR_Org               &  ctaOthe_NSR=ctaOthe_NSR+1l
                                                     ;;RTo_othe_R11=RTo_othe_R11+(reflectivZ/aRR)^(1/bRR)   &  ctaOthe_R11=ctaOthe_R11+1l   ;; SRB
                                                     
                                                     ;;rain_R11Full[tmp_col,tmp_row,3]=rain_R11Full[tmp_col,tmp_row,3]+(reflectivZ/aRR)^(1/bRR)  ;; SRB 
                                                     ;;nRai_R11Full[tmp_col,tmp_row,3]=nRai_R11Full[tmp_col,tmp_row,3]+1l                   ;; SRB
                                                  endif else begin
                                                     if raintype[col,fil] eq raintype_noRainValue then begin ;;No rain
                                                        RTo_noRa_NSR=RTo_noRa_NSR+0.   &   ctaNoRa_NSR=ctaNoRa_NSR+1l
                                                        ;;RTo_noRa_R11=RTo_noRa_R11+0.   &   ctaNoRa_R11=ctaNoRa_R11+1l                     ;; SRB
                                                        
                                                        ;;rain_R11Full[tmp_col,tmp_row,3]=rain_R11Full[tmp_col,tmp_row,3]+0.                ;; SRB
                                                        ;;nRai_R11Full[tmp_col,tmp_row,3]=nRai_R11Full[tmp_col,tmp_row,3]+1l                ;; SRB
                                                     endif
                                                  endelse
                                               endelse
                                            endelse
                                         endif    ;; end for missing values....
                                      endif       ;; if ctaC ne 0 and ctaR ne 0
                                   ;;endif       ;; end for pixels with height > 2.5km   ;; SRB
                                endif
                             endif
                          endfor
                       endif  ;;end for flag of accounting for a single FULL storm...

                       ;;print,'AT BOTTOM OF STORM CALCS: ss=',string(strtrim(ss,2)),' and ssST=',string(strtrim(ssST,2)), $
                       ;;      ' and area=',string(strtrim(area,2)),' and areaIsFull=',string(strtrim(areaIsFull,2))
                       ;;;;print,'     ctaStra_R11=',string(strtrim(ctaStra_R11,2)),' and ctaConv_R11=',string(strtrim(ctaConv_R11,2))
                       ;;print,'     ctaStra_NSR=',string(strtrim(ctaStra_NSR,2)),' and ctaConv_NSR=',string(strtrim(ctaConv_NSR,2))

                       ;;************************************************************************************************************   ;; SRB
                       ;;;;Romatschke and Houze 2011
                       ;;total_RainAll_R11=RTo_stra_R11+RTo_conv_R11+RTo_othe_R11+RTo_noRa_R11
                       ;;total_RainCSs_R11=RTo_stra_R11+RTo_conv_R11
                       ;;m_RainAll_R11=total_RainAll_R11/float(ctaStra_R11+ctaConv_R11+ctaOthe_R11+ctaNoRa_R11)  ;;mean rainfall for all pixels
                       ;;if ctaStra_R11 ne 0 then m_RainStrt_R11=RTo_stra_R11/ctaStra_R11 else m_RainStrt_R11=0. ;;mean stratiform rain
                       ;;if ctaConv_R11 ne 0 then m_RainConv_R11=RTo_conv_R11/ctaConv_R11 else m_RainConv_R11=0. ;;mean convective rain
                       ;;if ctaStra_R11 ne 0 or ctaConv_R11 ne 0 then $
                       ;;   m_RainCSs_R11=total_RainCSs_R11/float(ctaStra_R11+ctaConv_R11) else m_RainCSs_R11=-999.
                       ;;;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       ;;size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;vol_Rain_All_R11=total_RainAll_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Str_R11=RTo_stra_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Con_R11=RTo_conv_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       ;;tmp_statsRain_R11=[m_RainAll_R11,m_RainStrt_R11,m_RainConv_R11,vol_Rain_All_R11,vol_Rain_Str_R11,vol_Rain_Con_R11]

                       ;;************************************************************************************************************
                       ;;Near Surface Rain
                       total_RainAll_NSR=RTo_stra_NSR+RTo_conv_NSR+RTo_othe_NSR+RTo_noRa_NSR
                       total_RainCSs_NSR=RTo_stra_NSR+RTo_conv_NSR
                       m_RainAll_NSR=total_RainAll_NSR/float(ctaStra_NSR+ctaConv_NSR+ctaOthe_NSR+ctaNoRa_NSR)  ;;mean rainfall for all pixels
                       if ctaStra_NSR ne 0 then m_RainStrt_NSR=RTo_stra_NSR/ctaStra_NSR else m_RainStrt_NSR=0. ;;mean stratiform rain
                       if ctaConv_NSR ne 0 then m_RainConv_NSR=RTo_conv_NSR/ctaConv_NSR else m_RainConv_NSR=0. ;;mean convective rain
                       if ctaStra_NSR ne 0 or ctaConv_NSR ne 0 then $
                          m_RainCSs_NSR=total_RainCSs_NSR/float(ctaStra_NSR+ctaConv_NSR) else m_RainCSs_NSR=-999.
                       ;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       vol_Rain_All_NSR=total_RainAll_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Str_NSR=RTo_stra_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Con_NSR=RTo_conv_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       tmp_statsRain_NSR=[m_RainAll_NSR,m_RainStrt_NSR,m_RainConv_NSR,vol_Rain_All_NSR,vol_Rain_Str_NSR,vol_Rain_Con_NSR]

                       ;;*********************************************************************************************
                       ;;Now only for the Broad Strtiform part
                       ;;*******************************************************************************************
                       ;; ADDED _ST TO THESE VAR NAMES FOR CLARITY  ;; SRB 8/3/2018
                       ;; Z-R Method
                       ;;RTo_conv_R11_ST=0. & ctaConv_R11_ST=0l            ;;total convective rain   ;; SRB
                       ;;RTo_stra_R11_ST=0. & ctaStra_R11_ST=0l            ;;total stratiform rain   ;; SRB
                       ;;RTo_othe_R11_ST=0. & ctaOthe_R11_ST=0l            ;;total other rain        ;; SRB
                       ;;RTo_noRa_R11_ST=0. & ctaNoRa_R11_ST=0l            ;;total no_rain rain      ;; SRB

                      ;; ADDED _ST TO THESE VAR NAMES FOR CLARITY  ;; SRB 8/3/2018
                       ;; Near Surface Rain
                       RTo_conv_NSR_ST=0. & ctaConv_NSR_ST=0l ;;total convective rain
                       RTo_stra_NSR_ST=0. & ctaStra_NSR_ST=0l ;;total stratiform rain
                       RTo_othe_NSR_ST=0. & ctaOthe_NSR_ST=0l ;;total other rain
                       RTo_noRa_NSR_ST=0. & ctaNoRa_NSR_ST=0l ;;total no_rain rain

                       for i=0l,pixelsumST_ST-1l do begin 
                          col=dondeST_ST[i] mod nlonsFull 
                          fil=long(fix(dondeST_ST[i]/float(nlonsFull)))  

                          nearSrfR_Org=SrfRain_FULL[col,fil]

                          tmp_col=where(lonsC le lonsFull_sub[col],ctaC) ;;count the pixels
                          tmp_row=where(latsC le latsFull_sub[fil],ctaR) 

                          if ctaC ne 0 and ctaR ne 0 then begin 
                             colBS=(reverse(tmp_col))[0] 
                             rowBS=(reverse(tmp_row))[0] 
                             freq_Core[colBS,rowBS,3]=freq_Core[colBS,rowBS,3]+1l 

                             pila=reform(grid_storm_FULL[col,fil,*])  ;;New Corrected: Matching location within the full storm containing the core
                             w_CV=where(pila ne 0,ctaHgt)

                             if ctaHgt ge 1 then begin 
                                ;;distance between lowest pixel and ground
                                id_top1=(where(topo_lon ge lonsFull_sub[col]))[0]
                                id_top2=(where(topo_lat ge latsFull_sub[fil]))[0]

                                ;;%we sort out all the pixels that are higher than 2.5km above the ground
                                ;;if (hgts[w_CV[0]]-float(DEM[id_top1,id_top2])/1000.) le 2.5 then begin                              ;; SRB
                                   ;;locate the pixel in the nearest fine grid cell
                                   tmp_col=(where(float(lonsF) eq lonsFull_sub[col],ctaC))[0]
                                   tmp_row=(where(float(latsF) eq latsFull_sub[fil],ctaR))[0]

                                   if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                   ;;if ctaC eq 0 or ctaR eq 0 then stop
                                      
                                      ;;if nearSrfR_Org ne -9999.00 and refl_3D_FULL[col,fil,w_CV[0]] ne -999.0 then begin               ;; SRB
                                      if nearSrfR_Org ne SrfRain_fillValue then begin
                                         ;;reflectivZ=10^(refl_3D_FULL[col,fil,w_CV[0]]*0.1)  ;;%convert from dBZ to Z                   ;; SRB
                           
                                         rain_NSRCore[tmp_col,tmp_row,3]=rain_NSRCore[tmp_col,tmp_row,3]+nearSrfR_Org
                                         nRai_NSRCore[tmp_col,tmp_row,3]=nRai_NSRCore[tmp_col,tmp_row,3]+1l

                                         ;;here I create accumulated rain rate vectors
                                         if raintypeFULL[col,fil] eq CONV then begin ;;convective rain
                                            RTo_conv_NSR_ST=RTo_conv_NSR_ST+nearSrfR_Org                     &  ctaConv_NSR_ST=ctaConv_NSR_ST+1l
                                            ;;RTo_conv_R11_ST=RTo_conv_R11_ST+(reflectivZ/aRRc)^(1/bRRc)       &  ctaConv_R11_ST=ctaConv_R11_ST+1l   ;; SRB
 
                                            ;;rain_R11Core[tmp_col,tmp_row,3]=rain_R11Core[tmp_col,tmp_row,3]+(reflectivZ/aRRc)^(1/bRRc) ;; SRB
                                            ;;nRai_R11Core[tmp_col,tmp_row,3]=nRai_R11Core[tmp_col,tmp_row,3]+1l                         ;; SRB
                                         endif else begin
                                            if raintypeFULL[col,fil] eq STRA then begin ;;stratiform rain
                                               RTo_stra_NSR_ST=RTo_stra_NSR_ST+nearSrfR_Org                  &  ctaStra_NSR_ST=ctaStra_NSR_ST+1l
                                               ;;RTo_stra_R11_ST=RTo_stra_R11_ST+(reflectivZ/aRRs)^(1/bRRs)    &  ctaStra_R11_ST=ctaStra_R11_ST+1l   ;; SRB

                                               ;;rain_R11Core[tmp_col,tmp_row,3]=rain_R11Core[tmp_col,tmp_row,3]+(reflectivZ/aRRs)^(1/bRRs) ;; SRB
                                               ;;nRai_R11Core[tmp_col,tmp_row,3]=nRai_R11Core[tmp_col,tmp_row,3]+1l                      ;; SRB
                                            endif else begin
                                               if raintypeFULL[col,fil] ge OTHER then begin ;;other rain
                                                  RTo_othe_NSR_ST=RTo_othe_NSR_ST+nearSrfR_Org               &  ctaOthe_NSR_ST=ctaOthe_NSR_ST+1l
                                                  ;;RTo_othe_R11_ST=RTo_othe_R11_ST+(reflectivZ/aRR)^(1/bRR)   &  ctaOthe_R11_ST=ctaOthe_R11_ST+1l   ;; SRB

                                                  ;;rain_R11Core[tmp_col,tmp_row,3]=rain_R11Core[tmp_col,tmp_row,3]+(reflectivZ/aRR)^(1/bRR) ;; SRB
                                                  ;;nRai_R11Core[tmp_col,tmp_row,3]=nRai_R11Core[tmp_col,tmp_row,3]+1l                   ;; SRB
                                               endif else begin
                                                  if raintypeFULL[col,fil] eq raintype_noRainValue then begin ;;No rain
                                                     RTo_noRa_NSR_ST=RTo_noRa_NSR_ST+0.   &   ctaNoRa_NSR_ST=ctaNoRa_NSR_ST+1l
                                                     ;;RTo_noRa_R11_ST=RTo_noRa_R11_ST+0.   &   ctaNoRa_R11_ST=ctaNoRa_R11_ST+1l                     ;; SRB

                                                     ;;rain_R11Core[tmp_col,tmp_row,3]=rain_R11Core[tmp_col,tmp_row,3]+0.                ;; SRB
                                                     ;;nRai_R11Core[tmp_col,tmp_row,3]=nRai_R11Core[tmp_col,tmp_row,3]+1l                ;; SRB
                                                  endif
                                               endelse
                                            endelse
                                         endelse
                                      endif  ;; end for missing values....
                                   endif     ;; if ctaC ne 0 and ctaR ne 0
                                ;;endif     ;; end for pixels with height > 2.5km                                                     ;; SRB
                             endif
                          endif
                       endfor 

                       ;;************************************************************************************************************ ;; SRB
                       ;;Romatschke and Houze 2011
                       ;;total_RainAll_R11=RTo_stra_R11_ST+RTo_conv_R11_ST+RTo_othe_R11_ST+RTo_noRa_R11_ST
                       ;;total_RainCSs_R11=RTo_stra_R11_ST+RTo_conv_R11_ST
                       ;;m_RainAll_R11=total_RainAll_R11/float(ctaStra_R11_ST+ctaConv_R11_ST+ctaOthe_R11_ST+ctaNoRa_R11_ST)  ;;mean rainfall for all pixels
                       ;;if ctaStra_R11_ST ne 0 then m_RainStrt_R11=RTo_stra_R11_ST/ctaStra_R11_ST else m_RainStrt_R11=0. ;;mean stratiform rain
                       ;;if ctaConv_R11_ST ne 0 then m_RainConv_R11=RTo_conv_R11_ST/ctaConv_R11_ST else m_RainConv_R11=0. ;;mean convective rain
                       ;;if ctaStra_R11_ST ne 0 or ctaConv_R11_ST ne 0 then $
                       ;;   m_RainCSs_R11=total_RainCSs_R11/float(ctaStra_R11_ST+ctaConv_R11_ST) else m_RainCSs_R11=-999.
                       ;;;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       ;;size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;vol_Rain_All_R11=total_RainAll_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Str_R11=RTo_stra_R11_ST*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Con_R11=RTo_conv_R11_ST*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       ;;tmp_statsRain_R11_BS=[m_RainAll_R11,m_RainStrt_R11,m_RainConv_R11,vol_Rain_All_R11,vol_Rain_Str_R11,vol_Rain_Con_R11]

                       ;;************************************************************************************************************
                       ;;Near Surface Rain
                       total_RainAll_NSR=RTo_stra_NSR_ST+RTo_conv_NSR_ST+RTo_othe_NSR_ST+RTo_noRa_NSR_ST
                       total_RainCSs_NSR=RTo_stra_NSR_ST+RTo_conv_NSR_ST
                       m_RainAll_NSR=total_RainAll_NSR/float(ctaStra_NSR_ST+ctaConv_NSR_ST+ctaOthe_NSR_ST+ctaNoRa_NSR_ST)  ;;mean rainfall for all pixels
                       if ctaStra_NSR_ST ne 0 then m_RainStrt_NSR=RTo_stra_NSR_ST/ctaStra_NSR_ST else m_RainStrt_NSR=0. ;;mean stratiform rain
                       if ctaConv_NSR_ST ne 0 then m_RainConv_NSR=RTo_conv_NSR_ST/ctaConv_NSR_ST else m_RainConv_NSR=0. ;;mean convective rain
                       if ctaStra_NSR_ST ne 0 or ctaConv_NSR_ST ne 0 then $
                          m_RainCSs_NSR=total_RainCSs_NSR/float(ctaStra_NSR_ST+ctaConv_NSR_ST) else m_RainCSs_NSR=-999.
                       ;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       vol_Rain_All_NSR=total_RainAll_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Str_NSR=RTo_stra_NSR_ST*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Con_NSR=RTo_conv_NSR_ST*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       tmp_statsRain_NSR_BS=[m_RainAll_NSR,m_RainStrt_NSR,m_RainConv_NSR,vol_Rain_All_NSR,vol_Rain_Str_NSR,vol_Rain_Con_NSR]
                       
                       ;;*********************************************************************************************
                       ;;Here I calculate the CFAD count !!! for Full Storm!!! (only for one single Full storm)
                       if area ne areaIsFULL and cen_lon ne lonCIsFULL and cen_lat ne latCIsFull then begin
                          refl_SingleStorm=fltarr(nlonsFull,nlatsFull,nlevels)
                          refl_SingleStorm[*,*,*]=refl_3D_fillValue
                          refl_SingleStorm[w_idF]=refl_3D_FULL[w_idF]
                          if cta_Full ne npix_str[donde_BrdStr[ss]] then stop  ;; just to check! because this is in 3D!
           
                          ;;here count reflectivity for each pixel that compose the storm into a matrix of CFAD
                          for i=0l,cta_Full-1l do begin  
                             col_Refl=where(refl_CFAD eq round(refl_3D_FULL[w_idF[i]]),ctaZ) 
                             row_hgts=where(alts_CFAD eq hgts_3D_FULL[w_idF[i]],ctaH)  ;;here locate the height of the pixel
                             if ctaH ne 0 and ctaZ ne 0 then $
                                CFAD_Full[col_Refl,row_hgts,3]=CFAD_Full[col_Refl,row_hgts,3]+1l  else stop 
                          endfor
                       endif

                       ;;*********************************************************************************************
                       ;;Here I calculate the CFAD count !!! for Broad Stratiform component!!!
                       refl_SingleStorm=fltarr(nlonsFull,nlatsFull,nlevels)
                       refl_SingleStorm[*,*,*]=refl_3D_fillValue
                       refl_SingleStorm[w_idST]=refl_3D_FULL[w_idST]
                       if cta_ST ne npix_ST[donde_BrdStr2[ssST]] then stop  ;; just to check! because this is in 3D!

                       ;;here count reflectivity for each pixel that compose the storm into a matrix of CFAD
                       for i=0l,cta_ST-1 do begin  
                          col_Refl=where(refl_CFAD eq round(refl_3D_FULL[w_idST[i]]),ctaZ) 
                          row_hgts=where(alts_CFAD eq hgts_3D_FULL[w_idST[i]],ctaH)  ;;here locate the height of the pixel
                          if ctaH ne 0 and ctaZ ne 0 then $
                             CFAD_Core[col_Refl,row_hgts,3]=CFAD_Core[col_Refl,row_hgts,3]+1l  else stop 
                       endfor

                       areaIsFULL=area & lonCIsFULL=cen_lon & latCIsFull=cen_lat   ;;return area counter to avoid double count of fullStorm
                       
                       ;;************************************************************************************************************
                       ;;store the info for monthly_class directory
                       ;;************************************************************************************************************
                       ;; add storm number to output for correlation with mask in 3D volume
                       ;;info_BS=[info_BS,orbit+'.'+datetime+'.'+strtrim(string(s_idST),2)]
                       info_BS=[info_BS,orbit+'.'+datetime+'.'+strtrim(string(num_bsr),2)]

                       shape_Core_BS=[[shape_Core_BS],[tmp_shapeBS]]
                       shape_Full_BS=[[shape_Full_BS],[tmp_shape]]

                       rain_Core_BS=[[rain_Core_BS],[rain_momentBS]]
                       rain_Full_BS=[[rain_Full_BS],[rain_moment]]

                       rainTypeCore_BS=[[rainTypeCore_BS],[statsRain_BS]]
                       rainTypeFull_BS=[[rainTypeFull_BS],[statsRain]]

                       ;;store the info for the stat_class directory...;;************************************************************
                       rainCore_BS_NSR=[[rainCore_BS_NSR],[tmp_statsRain_NSR_BS]]
                       rainFull_BS_NSR=[[rainFull_BS_NSR],[tmp_statsRain_NSR]]
                       
                       ;;rainCore_BS_R11=[[rainCore_BS_R11],[tmp_statsRain_R11_BS]]            ;; SRB
                       ;;rainFull_BS_R11=[[rainFull_BS_R11],[tmp_statsRain_R11]]               ;; SRB

                       ;;************************************************************************************************************
                       undefine,dim_lonBS
                       undefine,dim_latBS
                       undefine,hgt_sumBS
                       undefine,grid_sum
                       undefine,donde
                       undefine,pixelsum
                       undefine,lon_sum
                       undefine,lat_sum
                       undefine,hgt_sum
                       undefine,tmp_shapeBS
                       undefine,tmp_shape
                       undefine,SrfRain_FULL
                       undefine,raintypeFULL
                       undefine,refl_3D_FULL
                       undefine,hgts_3D_FULL
                       undefine,grid_storm_FULL

                       undefine,rainBS
                       undefine,rain
                       undefine,rain_nomiss
                       undefine,stratconv
                       undefine,strats
                       undefine,convec
                       undefine,others
                       undefine,noRain
                       undefine,missin
                       undefine,pila

                       undefine,col_Refl
                       undefine,row_hgts
                       undefine,refl_SingleStorm
                       undefine,colBS
                       undefine,rowBS
                       undefine,tmp_col
                       undefine,tmp_row

                    endif      ;;endif area_BS gt thr_aST
                    undefine,s_idST
                    undefine,w_idST
                    undefine,singlestormgrid_ST
                    undefine,grid_sum_ST_ST
                    undefine,dondeST_ST
                    undefine,lonST
                    undefine,latST
                 endfor     ;;endfor loop through pixels that maybe are broad Stratiform
                 undefine,grid_ST
                 undefine,npix_ST
                 undefine,id_ST
                 undefine,searchNaN_ST
                 undefine,donde_BrdStr2
              endif      ;;endif for stratiforms areas (area_ST) > thr_aST
              undefine,total_lonFull
              undefine,d_lonsFull
              undefine,total_latFull
              undefine,d_latsFull
              undefine,singlestormgrid_Full
              undefine,lonsFull_sub
              undefine,latsFull_sub
              undefine,singlestormgridStratiform
              undefine,grid_sumST
              undefine,dondeST
              undefine,s_idF
              undefine,w_idF
           endfor        ;;endfor loop thru stratiform volumes greater than 1000 pixels
        endif            ;;end if cta_strt ne 0 and ctaBrdStr gt 0
        
        ;;*******************
        ;;*******************
        print,'Id-Convective'
        ;;*******************
        ;;*******************
        ;;here it is going to identify Convective subset
        where_convective=where(refl_3D ge thr_dbZ,cta_thr_dbZ) ;;check for all pixels with >=thr_dbZ convective
        delvar,where_convective

        ;;here I only choose convective pixels and from those only >=2 pixels with refl>=thr_dbZ 
        donde_Convec=where(npix_str ge 2l,ctaConvective) 

        ;;*********************************************************************************************
        ;;Identify Convective cores (DCC-WCC-DWC)
        ;; if cta_conv (num conv pixels in raintype) ne 0 and
        ;;    id_storm from first find_storm call not empty and
        ;;       number entries in refl_3D > thr_dbZ, >= 2
        if cta_conv ne 0 and id_storm[0] ne -999l and cta_thr_dbZ ge 2 and ctaConvective gt 0 then begin
           for ss=0l,ctaConvective-1 do begin ;;only goes thru the volumes having more than 25 pixels
              s_idF=id_storm[donde_Convec[ss]]
              w_idF=where(grid_storm eq s_idF,cta1)
              if cta1 ne npix_str[donde_Convec[ss]] then stop  ;; just to check!
              singlestormgrid=lonarr(nlons,nlats,nlevels)
              singlestormgrid[w_idF]=s_idF
              
              ;;Here I will subset the singleFULL storm found to reduce the size of the matrix
              total_lonFull=total(total(singlestormgrid,2),2) 
              d_lonsFull=where(total_lonFull gt 0l,nlonsFull) ;;this are the horizontal positions within a selcted contiguous area
              total_latFull=total(total(singlestormgrid,1),2) 
              d_latsFull=where(total_latFull gt 0l,nlatsFull)
              
              singlestormgrid_Full=lonarr(nlonsFull,nlatsFull,nlevels)
              singlestormgrid_Full=singlestormgrid[d_lonsFull,d_latsFull,*]
              w_idF=where(singlestormgrid_Full eq s_idF,cta_Full)
              lonsFull_sub=lons[d_lonsFull]
              latsFull_sub=lats[d_latsFull]
         
              undefine,singlestormgrid

              ;;Identify the convective pixels in storm (refl >= thr_dbZ and All convective pixels)
              singlestormgridConvective=singlestormgrid_Full
              singlestormgridConvective[where(refl_3D[d_lonsFull,d_latsFull,*] lt thr_dbZ)]=0l       ;;(refl >= thr_dbZ
              singlestormgridConvective[where(rain_type3D[d_lonsFull,d_latsFull,*] ne CONV)]=0l       ;;only convective pixels=20
               
              ;;2D horizontal array of number of Convective pixels within storm cluster
              grid_sumCV=total(singlestormgridConvective,3) 
              dondeCV=where(grid_sumCV gt 0,pixelsumCV) ;;pixelsum=number of pixels in a 2D proj

              if pixelsumCV ge 2 then begin ;;this is to avoid possible events with zero convective pixels
                 lonCV=total(total(singlestormgridConvective,2),2) 
                 lonC_CV=(max(lonsFull_sub[where(lonCV gt 0l)])+min(lonsFull_sub[where(lonCV gt 0l)]))/2. ;;center of Stratiform

                 latCV=total(total(singlestormgridConvective,1),2) 
                 latC_CV=(max(latsFull_sub[where(latCV gt 0l)])+min(latsFull_sub[where(latCV gt 0l)]))/2. ;;center of Stratiform

                 size_pixels=deg2km(pixDeg,lonC_CV,latC_CV)
                 area_CV=pixelsumCV*(size_pixels[0]*size_pixels[1]) ;;Convective area in km2  

                 hgt_sumCV=total(total(singlestormgridConvective,2),1)
                 dim_hgtCV=max(hgts[where(hgt_sumCV gt 0l)])-min(hgts[where(hgt_sumCV gt 0l)]) ;;vertical dimension in km
                 dim_topCV=max(hgts[where(hgt_sumCV gt 0l)])                                   ;;###############
                 dim_botCV=min(hgts[where(hgt_sumCV gt 0l)])                                   ;;###############
              endif else begin
                 area_CV=0.
                 dim_hgtCV=0.
                 dim_topCV=0.
              endelse

              ;;check if there are pixels within the storm with area>=thr_aCV or height>=thr_hCV
              if area_CV ge thr_aCV or dim_topCV ge thr_hCV then begin ;;original Ulli's classification
                 ;;Does the convective area belongs to a single storm? (acomplish contiguous pixel condition)
                 ;;Again, uses the convective pixels to identify contiguous pixels         
                 findStormNew,refl_3D=singlestormgridConvective,$
                              ;;refl_3D_fillValue=refl_3D_fillValue,$
                              id_storm=id_CV,$
                              npix_str=npix_CV,$
                              grid_storm=grid_CV
                 searchNaN_CV=where(grid_CV lt 0, nanCnt)
                 if nanCnt gt 0 then grid_CV[searchNaN_CV]=0l ;;%set NaNs to zero in the  grid matrix 
                 
                 ;;here I only choose pixel areas with more than >2 pixles (to avoid smaller areas)
                 donde_Convec2=where(npix_CV ge 2l,ctaConvective2) 

                 areaIsFULL=0. & lonCIsFULL=0. & latCIsFull=0. ;;this is to locate only ONE full storm
                 
                 for ssCV=0l,ctaConvective2-1 do begin ;;only goes thru volumes greater than 25 pixels (speed the process)
                    s_idCV=id_CV[donde_Convec2[ssCV]]
                    w_idCV=where(grid_CV eq s_idCV,cta_CV)
                    if cta_CV ne npix_CV[donde_Convec2[ssCV]] then stop  ;; just to check!

                    singlestormgrid_CV=lonarr(nlonsFull,nlatsFull,nlevels)
                    singlestormgrid_CV[w_idCV]=s_idCV

                    ;;2D horizontal array of number of Stratiform pixels within storm
                    grid_sum_CV_CV=total(singlestormgrid_CV,3) 
                    dondeCV_CV=where(grid_sum_CV_CV gt 0,pixelsumCV_CV) ;;pixelsum=number of pixels in a 2D proj

                    lonCV=total(total(singlestormgrid_CV,2),2) 
                    lonC_CV=(max(lonsFull_sub[where(lonCV gt 0l)])+min(lonsFull_sub[where(lonCV gt 0l)]))/2. ;;center of Convective

                    latCV=total(total(singlestormgrid_CV,1),2) 
                    latC_CV=(max(latsFull_sub[where(latCV gt 0l)])+min(latsFull_sub[where(latCV gt 0l)]))/2. ;;center of Convective

                    size_pixels=deg2km(pixDeg,lonC_CV,latC_CV)
                    area_CV=pixelsumCV_CV*(size_pixels[0]*size_pixels[1]) ;;Convective area in km  

                    hgt_sumCV=total(total(singlestormgrid_CV,2),1)
                    dim_hgtCV=max(hgts[where(hgt_sumCV gt 0l)])-min(hgts[where(hgt_sumCV gt 0l)]) ;;vertical dimension km
                    dim_topCV=max(hgts[where(hgt_sumCV gt 0l)])                                   ;;###############
                    dim_botCV=min(hgts[where(hgt_sumCV gt 0l)])                                   ;;###############
                    
                    ;; ADDED TEST FOR DCC TO BE BASED ON dim_topCV AND dim_hgtCV (old way just checked dim_topCV)
                    ;;now identify the real contiguous areas belonging to a Deep or Wide convective event
                    if pixelsumCV_CV ge 2 and (area_CV ge thr_aCV or (dim_topCV ge thr_hCV and dim_hgtCV ge thr_dCV) ) then begin ;;orig Ulli's classification

                       ;; add core to appropriate mask
                       if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective
                          num_dwc = num_dwc + 1
                          ;; add core id to dwc_mask if makeCoreMasks is set
                          if makeCoreMasks then begin
                             dwc_mask_sub = dwc_mask[d_lonsFull,d_latsFull,0]
                             dwc_mask_sub[dondeCV_CV]=num_dwc
                             dwc_mask[d_lonsFull,d_latsFull,0] = dwc_mask_sub
                             undefine,dwc_mask_sub
                          endif 
                       endif else begin
                          if area_CV ge thr_aCV then begin
                             num_wcc = num_wcc + 1
                             ;; add core id to wcc_mask if makeCoreMasks is set
                             if makeCoreMasks then begin
                                wcc_mask_sub = wcc_mask[d_lonsFull,d_latsFull,0]
                                wcc_mask_sub[dondeCV_CV]=num_wcc
                                wcc_mask[d_lonsFull,d_latsFull,0] = wcc_mask_sub
                                undefine,wcc_mask_sub
                             endif 
                          endif else begin
                             if dim_topCV ge thr_hCV then begin
                                num_dcc = num_dcc + 1
                                ;; add core id to dcc_mask if makeCoreMasks is set
                                if makeCoreMasks then begin
                                   dcc_mask_sub = dcc_mask[d_lonsFull,d_latsFull,0]
                                   dcc_mask_sub[dondeCV_CV]=num_dcc
                                   dcc_mask[d_lonsFull,d_latsFull,0] = dcc_mask_sub
                                   undefine,dcc_mask_sub
                                endif 
                             endif 
                          endelse
                       endelse

                       ;;calculate the dimensions for Convective pixels in the cluster
                       ;;*************************************************************************************
                       dim_lonCV=(max(lonsFull_sub[where(lonCV gt 0l)])-min(lonsFull_sub[where(lonCV gt 0l)]))+pixDeg
                       dim_latCV=(max(latsFull_sub[where(latCV gt 0l)])-min(latsFull_sub[where(latCV gt 0l)]))+pixDeg

                       ;;calculates the elevation of the terrain for the center of the storm  !21
                       id_top1=(where(topo_lon ge lonC_CV))[0] & id_top2=(where(topo_lat ge latC_CV))[0]
                       terr_hgtCV=DEM[id_top1,id_top2]                            ;;elevation in meters
                       if terr_hgtCV eq 0 then land_oceanCV=0 else land_oceanCV=1 ;;ocean=0 or land=1!

                       tmp_shapeCV=[lonC_CV,latC_CV,area_CV,dim_topCV,dim_botCV,dim_lonCV,dim_latCV,terr_hgtCV,land_oceanCV]
                       
                       ;;now it calculate the dimensions for full storm (all pixels in the cluster)
                       ;;******************************************************************************************
                       grid_sum=total(singlestormgrid_Full,3) ;;2D horiz array of number of all pixels within storm
                       donde=where(grid_sum gt 0,pixelsum)    ;;pixelsum=number of pixels in a 2D proj

                       ;; increment num_storm
                       num_storm = num_storm + 1

                       ;; add storm to storm_mask if makeCoreMasks is set
                       if makeCoreMasks then begin
                          storm_mask_sub = storm_mask[d_lonsFull,d_latsFull,0]
                          storm_mask_sub[donde]=num_storm
                          storm_mask[d_lonsFull,d_latsFull,0] = storm_mask_sub
                          undefine,storm_mask_sub
                       endif 

                       lon_sum=total(total(singlestormgrid_Full,2),2)
                       cen_lon=(max(lonsFull_sub[where(lon_sum gt 0l)])+min(lonsFull_sub[where(lon_sum gt 0l)]))/2.     ;;longitude coordinate
                       dim_lon=(max(lonsFull_sub[where(lon_sum gt 0l)])-min(lonsFull_sub[where(lon_sum gt 0l)]))+pixDeg ;;longitudinal extents

                       lat_sum=total(total(singlestormgrid_Full,1),2)  
                       cen_lat=(max(latsFull_sub[where(lat_sum gt 0l)])+min(latsFull_sub[where(lat_sum gt 0l)]))/2.     ;;latitude coordinate
                       dim_lat=(max(latsFull_sub[where(lat_sum gt 0l)])-min(latsFull_sub[where(lat_sum gt 0l)]))+pixDeg ;;latitudinal extents

                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       area=pixelsum*size_pixels[0]*size_pixels[1] ;;in km horizontal area of selected storm
     
                       hgt_sum=total(total(singlestormgrid_Full,2),1)
                       dim_hgt=(max(hgts[where(hgt_sum gt 0l)])-min(hgts[where(hgt_sum gt 0l)])) ;;%vertical extent in km
                       dim_top=max(hgts[where(hgt_sum gt 0l)])                                   ;;###############
                       dim_bot=min(hgts[where(hgt_sum gt 0l)])                                   ;;###############

                       ;;%calculates the elevation of the terrain for the center of the storm  !21
                       id_top1=(where(topo_lon ge cen_lon))[0] & id_top2=(where(topo_lat ge cen_lat))[0]
                       terr_hgt=DEM[id_top1,id_top2]                        ;;elevation in meters
                       if terr_hgt eq 0 then land_ocean=0 else land_ocean=1 ;;mask  ocean=0 or land=1!
                       
                       tmp_shape=[cen_lon,cen_lat,area,dim_top,dim_bot,dim_lon,dim_lat,terr_hgt,land_ocean]

                       ;;******************************************************************************************
                       ;;Statistics for rain types and rain rates *************************************************
                       ;;Statistics over Convective subset!
                       SrfRain_FULL=SrfRain[d_lonsFull,d_latsFull,*]
                       raintypeFULL=raintype[d_lonsFull,d_latsFull,*]
                       refl_3D_FULL=refl_3D[d_lonsFull,d_latsFull,*]
                       hgts_3D_FULL=hgts_3D[d_lonsFull,d_latsFull,*]

                       stratconv=raintypeFULL[dondeCV_CV]                     ;;type of rain in each 2D pixel that compose the storm
                       strats=where(stratconv eq STRA,ctaStr)                 ;;stratiform
                       convec=where(stratconv eq CONV,ctaCon)                 ;;convective
                       others=where(stratconv ge OTHER,ctaOth)                ;;other type
                       noRain=where(stratconv eq raintype_noRainValue,ctaNoR) ;;no rain
                       missin=where(stratconv eq raintype_fillValue,ctaMis)   ;;missing value
                       
                       rainCV=SrfRain_FULL[dondeCV_CV] ;;this is based on Near Surface Rain
                       rain_nomiss=where(rainCV ne SrfRain_fillValue,Rmiss)

                       ;;here I calculate moments of simple rain within storm
                       if Rmiss ge 2 then rain_momentCV=[mean(rainCV[rain_nomiss]),stdev(rainCV[rain_nomiss]),$
                                                         max(rainCV[rain_nomiss]),min(rainCV[rain_nomiss]),float(Rmiss),$
                                                         float(ctaStr),float(ctaCon)] $
                       else rain_momentCV=[-9999.,-9999.,-9999.,-9999.,-9999.,-9999.,-9999.]

                       ;;here I calculate rainrate sums in [mm/hr]
                       if ctaStr ne 0 then RTo_stra=total(rainCV[strats]) else RTo_stra=0. ;;total stratiform rain
                       if ctaCon ne 0 then RTo_conv=total(rainCV[convec]) else RTo_conv=0. ;;total convective rain
                       if ctaOth ne 0 then RTo_othe=total(rainCV[others]) else RTo_othe=0. ;;total other rain
                       if ctaNoR ne 0 then RTo_noRa=total(rainCV[noRain]) else RTo_noRa=0. ;;total no_rain rain
                       total_RainAll=RTo_stra+RTo_conv+RTo_othe+RTo_noRa
                       total_RainCSs=RTo_stra+RTo_conv

                       m_RainAll=total_RainAll/(ctaStr+ctaCon+ctaOth+ctaNoR+ctaMis)      ;;mean rainfall all pixels
                       if ctaStr ne 0 then m_RainStrt=RTo_stra/ctaStr else m_RainStrt=0. ;;mean stratiform rain
                       if ctaCon ne 0 then m_RainConv=RTo_conv/ctaCon else m_RainConv=0. ;;mean convective rain
                       if ctaStr ne 0 or ctaCon ne 0 then m_RainCSs=total_RainCSs/(ctaStr+ctaCon) $
                       else m_RainCSs=-999.

                       ;;here calculate volumetric values of rain  ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,lonC_CV,latC_CV)
                       vol_Rain_All=total_RainAll*(size_pixels[0]*size_pixels[1])/secsPerHr  
                       vol_Rain_Str=RTo_stra*(size_pixels[0]*size_pixels[1])/secsPerHr 
                       vol_Rain_Con=RTo_conv*(size_pixels[0]*size_pixels[1])/secsPerHr
   
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       statsRain_CV=[m_RainAll,m_RainStrt,m_RainConv,vol_Rain_All,vol_Rain_Str,vol_Rain_Con]
 
                       ;;Statistics for rain types and rain rates **********************************************
                       ;;***************************************************************************************
                       stratconv=raintypeFULL[donde]                          ;;type of rain in each 2D pixel that compose the storm
                       strats=where(stratconv eq STRA,ctaStr)                 ;;stratiform
                       convec=where(stratconv eq CONV,ctaCon)                 ;;convective
                       others=where(stratconv ge OTHER,ctaOth)                ;;other type
                       noRain=where(stratconv eq raintype_noRainValue,ctaNoR) ;;no rain
                       missin=where(stratconv eq raintype_fillValue,ctaMis)   ;;missing value

                       rain=SrfRain_FULL[donde] ;;this is based on Near Surface Rain
                       rain_nomiss=where(rain ne SrfRain_fillValue,Rmiss)
                       
                       ;;here I calculate moments of simple rain within storm
                       if Rmiss ge 2 then rain_moment=[mean(rain[rain_nomiss]),stdev(rain[rain_nomiss]),$
                                                       max(rain[rain_nomiss]),min(rain[rain_nomiss]),float(Rmiss),$
                                                       float(ctaStr),float(ctaCon)] $
                       else rain_moment=[-9999.,-9999.,-9999.,-9999.,-9999.,-9999.,-9999.]

                       ;;here I calculate rainrate sums in [mm/hr]
                       if ctaStr ne 0 then RTo_stra=total(rain[strats]) else RTo_stra=0. ;;total stratiform rain
                       if ctaCon ne 0 then RTo_conv=total(rain[convec]) else RTo_conv=0. ;;total convective rain
                       if ctaOth ne 0 then RTo_othe=total(rain[others]) else RTo_othe=0. ;;total other rain
                       if ctaNoR ne 0 then RTo_noRa=total(rain[noRain]) else RTo_noRa=0. ;;total no_rain rain
                       total_RainAll=RTo_stra+RTo_conv+RTo_othe+RTo_noRa
                       total_RainCSs=RTo_stra+RTo_conv

                       m_RainAll=total_RainAll/(ctaStr+ctaCon+ctaOth+ctaNoR+ctaMis)      ;;mean rainfall all pixels
                       if ctaStr ne 0 then m_RainStrt=RTo_stra/ctaStr else m_RainStrt=0. ;;mean stratiform rain
                       if ctaCon ne 0 then m_RainConv=RTo_conv/ctaCon else m_RainConv=0. ;;mean convective rain
                       if ctaStr ne 0 or ctaCon ne 0 then m_RainCSs=total_RainCSs/(ctaStr+ctaCon) $
                       else m_RainCSs=-999.

                       ;;here calculate volumetric values of rain  ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;area=pixelsum*size_pixels[0]*size_pixels[1] ;;in km horizontal area of selected storm

                       vol_Rain_All=total_RainAll*(size_pixels[0]*size_pixels[1])/secsPerHr  
                       vol_Rain_Str=RTo_stra*(size_pixels[0]*size_pixels[1])/secsPerHr 
                       vol_Rain_Con=RTo_conv*(size_pixels[0]*size_pixels[1])/secsPerHr

                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       statsRain=[m_RainAll,m_RainStrt,m_RainConv,vol_Rain_All,vol_Rain_Str,vol_Rain_Con]

                       ;;**********************************************************************************************
                       grid_storm_FULL=lonarr(nlonsFull,nlatsFull,nlevels)
                       grid_storm_FULL=grid_storm[d_lonsFull,d_latsFull,*]

                       ;;----------------------------------------------------------------------------------------------
                       ;;MOVED THIS CODE INTO IF LOOP WHERE IT NEEDS TO BE TO GET RID OF NaN OUTPUTS
                       ;;Statistics for rain rates Calculated using Z-R method 
                       ;;we need to get all the needed information for the lowest pixels in each storm
                       ;;locate the indices of the lowest pixels with data   (Full Storm) (Romatschke and Houze 2011
                       ;;;;RTo_conv_R11=0. & ctaConv_R11=0l ;;total convective rain       ;; SRB
                       ;;;;RTo_stra_R11=0. & ctaStra_R11=0l ;;total stratiform rain       ;; SRB
                       ;;;;RTo_othe_R11=0. & ctaOthe_R11=0l ;;total other rain            ;; SRB
                       ;;;;RTo_noRa_R11=0. & ctaNoRa_R11=0l ;;total no_rain rain          ;; SRB
                
                       ;;Near Surface Rain
                       ;;RTo_conv_NSR=0. & ctaConv_NSR=0l ;;total convective rain
                       ;;RTo_stra_NSR=0. & ctaStra_NSR=0l ;;total stratiform rain
                       ;;RTo_othe_NSR=0. & ctaOthe_NSR=0l ;;total other rain
                       ;;RTo_noRa_NSR=0. & ctaNoRa_NSR=0l ;;total no_rain rain
                       ;;----------------------------------------------------------------------------------------------

                       ;;print,'BEFORE STORM CALCS: ss=',string(strtrim(ss,2)),' and ssCV=',string(strtrim(ssCV,2)), $
                       ;;      ' and area=',string(strtrim(area,2)),' and areaIsFULL=',string(strtrim(areaIsFULL,2))

                       ;;here also is accounted n_pixels in matrix matching with storm location and locate indices (lowest pixels           
                       ;;this is only evaluated ONCE! (accouting for a single FULL storm)
                       if area ne areaIsFULL and cen_lon ne lonCIsFULL and cen_lat ne latCIsFull then begin

                          ;;print,'     RUNNING THROUGH CALCS FOR STORM - SHOULD BE ONLY ONCE'

                          ;;Statistics for rain rates Calculated using Z-R method 
                          ;;RTo_conv_R11=0. & ctaConv_R11=0l ;;total convective rain       ;; SRB
                          ;;RTo_stra_R11=0. & ctaStra_R11=0l ;;total stratiform rain       ;; SRB
                          ;;RTo_othe_R11=0. & ctaOthe_R11=0l ;;total other rain            ;; SRB
                          ;;RTo_noRa_R11=0. & ctaNoRa_R11=0l ;;total no_rain rain          ;; SRB

                          ;;Statistics for rain rates Calculated using Near Surface Rain
                          RTo_conv_NSR=0. & ctaConv_NSR=0l ;;total convective rain
                          RTo_stra_NSR=0. & ctaStra_NSR=0l ;;total stratiform rain
                          RTo_othe_NSR=0. & ctaOthe_NSR=0l ;;total other rain
                          RTo_noRa_NSR=0. & ctaNoRa_NSR=0l ;;total no_rain rain

                          for i=0l,pixelsum-1l do begin   
                             col=donde[i] mod nlonsFull               ;;column ID of the pixel
                             fil=long(fix(donde[i]/float(nlonsFull))) ;;row    ID of the pixel
                      
                             nearSrfR_Org=SrfRain_FULL[col,fil]
                             
                             t_col=where(lonsC le lonsFull_sub[col],ctaC) ;;count the pixels
                             t_row=where(latsC le latsFull_sub[fil],ctaR) 
                             if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                colCV=(reverse(t_col))[0] 
                                rowCV=(reverse(t_row))[0] 
                                if area_CV ge thr_aCV and dim_topCV ge thr_hCV then $
                                   freq_Full[colCV,rowCV,2]=freq_Full[colCV,rowCV,2]+1l else $
                                      if area_CV ge thr_aCV then freq_Full[colCV,rowCV,1]=freq_Full[colCV,rowCV,1]+1l $
                                      else if dim_topCV ge thr_hCV then freq_Full[colCV,rowCV,0]=freq_Full[colCV,rowCV,0]+1l 

                                ;;locate indices (x,y,z location of individual pixels within the FULL storm)
                                pila=reform(grid_storm_FULL[col,fil,*])
                                w_CV=where(pila ne 0,ctaHgt)
                                if ctaHgt ge 1 then begin 
                                   ;;distance between lowest pixel and ground
                                   id_top1=(where(topo_lon ge lonsFull_sub[col]))[0]
                                   id_top2=(where(topo_lat ge latsFull_sub[fil]))[0]

                                   ;;%we sort out all the pixels that are higher than 2.5km above the ground
                                   ;;if (hgts[w_CV[0]]-float(DEM[id_top1,id_top2])/1000.) le 2.5 then begin                     ;; SRB
                                      ;;locate the pixel in the nearest fine grid cell
                                      tmp_col=(where(float(lonsF) eq lonsFull_sub[col],ctaC))[0]
                                      tmp_row=(where(float(latsF) eq latsFull_sub[fil],ctaR))[0]
                               
                                      if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                      ;;if ctaC eq 0 or ctaR eq 0 then stop
                               
                                         ;;if nearSrfR_Org ne -9999.00 and refl_3D_FULL[col,fil,w_CV[0]] ne -999.0 then begin      ;; SRB
                                         if nearSrfR_Org ne SrfRain_fillValue then begin
                                            ;;reflectivZ=10^(refl_3D_FULL[col,fil,w_CV[0]]*0.1)  ;;%convert from dBZ to Z          ;; SRB
                                  
                                            ;;here I create accumulated matrices with rain rate for each type of convective system
                                            if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective
                                               rain_NSRFull[tmp_col,tmp_row,2]=rain_NSRFull[tmp_col,tmp_row,2]+nearSrfR_Org
                                               nRai_NSRFull[tmp_col,tmp_row,2]=nRai_NSRFull[tmp_col,tmp_row,2]+1l
                                            endif else begin
                                               if area_CV ge thr_aCV then begin
                                                  rain_NSRFull[tmp_col,tmp_row,1]=rain_NSRFull[tmp_col,tmp_row,1]+nearSrfR_Org
                                                  nRai_NSRFull[tmp_col,tmp_row,1]=nRai_NSRFull[tmp_col,tmp_row,1]+1l
                                               endif else begin
                                                  if dim_topCV ge thr_hCV then begin
                                                     rain_NSRFull[tmp_col,tmp_row,0]=rain_NSRFull[tmp_col,tmp_row,0]+nearSrfR_Org
                                                     nRai_NSRFull[tmp_col,tmp_row,0]=nRai_NSRFull[tmp_col,tmp_row,0]+1l
                                                  endif 
                                               endelse
                                            endelse
                                            
                                            ;;here I create accumulated rain rate vectors
                                            if raintypeFULL[col,fil] eq CONV then begin ;;convective rain
                                               RTo_conv_NSR=RTo_conv_NSR+nearSrfR_Org                     &  ctaConv_NSR=ctaConv_NSR+1l
                                               ;;RTo_conv_R11=RTo_conv_R11+(reflectivZ/aRRc)^(1/bRRc)       &  ctaConv_R11=ctaConv_R11+1l ;; SRB

                                               ;;tmp_rain=(reflectivZ/aRRc)^(1/bRRc)                                                      ;; SRB
                                               ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective     ;; SRB
                                               ;;   rain_R11Full[tmp_col,tmp_row,2]=rain_R11Full[tmp_col,tmp_row,2]+tmp_rain              ;; SRB
                                               ;;   nRai_R11Full[tmp_col,tmp_row,2]=nRai_R11Full[tmp_col,tmp_row,2]+1l                    ;; SRB
                                               ;;endif else begin                                                                         ;; SRB
                                               ;;   if area_CV ge thr_aCV then begin                                                      ;; SRB
                                               ;;      rain_R11Full[tmp_col,tmp_row,1]=rain_R11Full[tmp_col,tmp_row,1]+tmp_rain           ;; SRB
                                               ;;      nRai_R11Full[tmp_col,tmp_row,1]=nRai_R11Full[tmp_col,tmp_row,1]+1l                 ;; SRB
                                               ;;   endif else begin                                                                      ;; SRB
                                               ;;      if dim_topCV ge thr_hCV then begin                                                 ;; SRB
                                               ;;         rain_R11Full[tmp_col,tmp_row,0]=rain_R11Full[tmp_col,tmp_row,0]+tmp_rain        ;; SRB
                                               ;;         nRai_R11Full[tmp_col,tmp_row,0]=nRai_R11Full[tmp_col,tmp_row,0]+1l              ;; SRB
                                               ;;      endif                                                                              ;; SRB
                                               ;;   endelse                                                                               ;; SRB
                                               ;;endelse                                                                                  ;; SRB

                                            endif else begin
                                               if raintypeFULL[col,fil] eq STRA then begin ;;stratiform rain
                                                  RTo_stra_NSR=RTo_stra_NSR+nearSrfR_Org                  &  ctaStra_NSR=ctaStra_NSR+1l
                                                  ;;RTo_stra_R11=RTo_stra_R11+(reflectivZ/aRRs)^(1/bRRs)    &  ctaStra_R11=ctaStra_R11+1l ;; SRB

                                                  ;;tmp_rain=(reflectivZ/aRRs)^(1/bRRs)                                                   ;; SRB
                                                  ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective  ;; SRB
                                                  ;;   rain_R11Full[tmp_col,tmp_row,2]=rain_R11Full[tmp_col,tmp_row,2]+tmp_rain           ;; SRB
                                                  ;;   nRai_R11Full[tmp_col,tmp_row,2]=nRai_R11Full[tmp_col,tmp_row,2]+1l                 ;; SRB
                                                  ;;endif else begin                                                                      ;; SRB
                                                  ;;   if area_CV ge thr_aCV then begin                                                   ;; SRB
                                                  ;;      rain_R11Full[tmp_col,tmp_row,1]=rain_R11Full[tmp_col,tmp_row,1]+tmp_rain        ;; SRB
                                                  ;;      nRai_R11Full[tmp_col,tmp_row,1]=nRai_R11Full[tmp_col,tmp_row,1]+1l              ;; SRB
                                                  ;;   endif else begin                                                                   ;; SRB
                                                  ;;      if dim_topCV ge thr_hCV then begin                                              ;; SRB
                                                  ;;         rain_R11Full[tmp_col,tmp_row,0]=rain_R11Full[tmp_col,tmp_row,0]+tmp_rain     ;; SRB
                                                  ;;         nRai_R11Full[tmp_col,tmp_row,0]=nRai_R11Full[tmp_col,tmp_row,0]+1l           ;; SRB
                                                  ;;      endif                                                                           ;; SRB
                                                  ;;   endelse                                                                            ;; SRB
                                                  ;;endelse                                                                               ;; SRB

                                               endif else begin
                                                  if raintypeFULL[col,fil] ge OTHER then begin ;;other rain
                                                     RTo_othe_NSR=RTo_othe_NSR+nearSrfR_Org               &  ctaOthe_NSR=ctaOthe_NSR+1l
                                                     ;;RTo_othe_R11=RTo_othe_R11+(reflectivZ/aRR)^(1/bRR)   &  ctaOthe_R11=ctaOthe_R11+1l    ;; SRB
                                                     
                                                     ;;tmp_rain=(reflectivZ/aRR)^(1/bRR)                                                     ;; SRB
                                                     ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective  ;; SRB
                                                     ;;   rain_R11Full[tmp_col,tmp_row,2]=rain_R11Full[tmp_col,tmp_row,2]+tmp_rain           ;; SRB
                                                     ;;   nRai_R11Full[tmp_col,tmp_row,2]=nRai_R11Full[tmp_col,tmp_row,2]+1l                 ;; SRB
                                                     ;;endif else begin                                                                      ;; SRB
                                                     ;;   if area_CV ge thr_aCV then begin                                                   ;; SRB
                                                     ;;      rain_R11Full[tmp_col,tmp_row,1]=rain_R11Full[tmp_col,tmp_row,1]+tmp_rain        ;; SRB
                                                     ;;      nRai_R11Full[tmp_col,tmp_row,1]=nRai_R11Full[tmp_col,tmp_row,1]+1l              ;; SRB
                                                     ;;   endif else begin                                                                   ;; SRB
                                                     ;;      if dim_topCV ge thr_hCV then begin                                              ;; SRB
                                                     ;;         rain_R11Full[tmp_col,tmp_row,0]=rain_R11Full[tmp_col,tmp_row,0]+tmp_rain     ;; SRB
                                                     ;;         nRai_R11Full[tmp_col,tmp_row,0]=nRai_R11Full[tmp_col,tmp_row,0]+1l           ;; SRB
                                                     ;;      endif                                                                           ;; SRB
                                                     ;;   endelse                                                                            ;; SRB
                                                     ;;endelse                                                                               ;; SRB

                                                  endif else begin
                                                     if raintype[col,fil] eq raintype_noRainValue then begin ;;No rain
                                                        RTo_noRa_NSR=RTo_noRa_NSR+0.   &   ctaNoRa_NSR=ctaNoRa_NSR+1l
                                                        ;;RTo_noRa_R11=RTo_noRa_R11+0.   &   ctaNoRa_R11=ctaNoRa_R11+1l                        ;; SRB
                                                        
                                                        ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective ;; SRB
                                                        ;;   rain_R11Full[tmp_col,tmp_row,2]=rain_R11Full[tmp_col,tmp_row,2]+0                 ;; SRB
                                                        ;;   nRai_R11Full[tmp_col,tmp_row,2]=nRai_R11Full[tmp_col,tmp_row,2]+1l                ;; SRB
                                                        ;;endif else begin                                                                     ;; SRB
                                                        ;;   if area_CV ge thr_aCV then begin                                                  ;; SRB
                                                        ;;      rain_R11Full[tmp_col,tmp_row,1]=rain_R11Full[tmp_col,tmp_row,1]+0              ;; SRB
                                                        ;;      nRai_R11Full[tmp_col,tmp_row,1]=nRai_R11Full[tmp_col,tmp_row,1]+1l             ;; SRB
                                                        ;;   endif else begin                                                                  ;; SRB
                                                        ;;      if dim_topCV ge thr_hCV then begin                                             ;; SRB
                                                        ;;         rain_R11Full[tmp_col,tmp_row,0]=rain_R11Full[tmp_col,tmp_row,0]+0           ;; SRB
                                                        ;;         nRai_R11Full[tmp_col,tmp_row,0]=nRai_R11Full[tmp_col,tmp_row,0]+1l          ;; SRB
                                                        ;;      endif                                                                          ;; SRB
                                                        ;;   endelse                                                                           ;; SRB
                                                        ;;endelse                                                                              ;; SRB
                                                     endif
                                                  endelse
                                               endelse
                                            endelse
                                         endif    ;; end for missing values....
                                      endif       ;; if ctaC ne 0 and ctaR ne 0
                                   ;;endif       ;; end for pixels with height > 2.5km
                                endif
                             endif
                          endfor
                       endif  ;;end for flag of accounting for a single FULL storm...

                       ;;print,'AT BOTTOM OF STORM CALCS: ss=',string(strtrim(ss,2)),' and ssCV=',string(strtrim(ssCV,2)), $
                       ;;      ' and area=',string(strtrim(area,2)),' and areaIsFull=',string(strtrim(areaIsFull,2))
                       ;;;;print,'     ctaStra_R11=',string(strtrim(ctaStra_R11,2)),' and ctaConv_R11=',string(strtrim(ctaConv_R11,2))         ;; SRB
                       ;;print,'     ctaStra_NSR=',string(strtrim(ctaStra_NSR,2)),' and ctaConv_NSR=',string(strtrim(ctaConv_NSR,2))
                       
                       ;;************************************************************************************************************        ;; SRB
                       ;;;;Romatschke and Houze 2011
                       ;;total_RainAll_R11=RTo_stra_R11+RTo_conv_R11+RTo_othe_R11+RTo_noRa_R11
                       ;;total_RainCSs_R11=RTo_stra_R11+RTo_conv_R11
                       ;;m_RainAll_R11=total_RainAll_R11/float(ctaStra_R11+ctaConv_R11+ctaOthe_R11+ctaNoRa_R11)  ;;mean rainfall for all pixels
                       ;;if ctaStra_R11 ne 0 then m_RainStrt_R11=RTo_stra_R11/ctaStra_R11 else m_RainStrt_R11=0. ;;mean stratiform rain
                       ;;if ctaConv_R11 ne 0 then m_RainConv_R11=RTo_conv_R11/ctaConv_R11 else m_RainConv_R11=0. ;;mean convective rain
                       ;;if ctaStra_R11 ne 0 or ctaConv_R11 ne 0 then $
                       ;;   m_RainCSs_R11=total_RainCSs_R11/float(ctaStra_R11+ctaConv_R11) else m_RainCSs_R11=-999.
                       ;;;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       ;;size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;vol_Rain_All_R11=total_RainAll_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Str_R11=RTo_stra_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Con_R11=RTo_conv_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       ;;tmp_statsRain_R11=[m_RainAll_R11,m_RainStrt_R11,m_RainConv_R11,vol_Rain_All_R11,vol_Rain_Str_R11,vol_Rain_Con_R11]

                       ;;************************************************************************************************************
                       ;;Near Surface Rain
                       total_RainAll_NSR=RTo_stra_NSR+RTo_conv_NSR+RTo_othe_NSR+RTo_noRa_NSR
                       total_RainCSs_NSR=RTo_stra_NSR+RTo_conv_NSR
                       m_RainAll_NSR=total_RainAll_NSR/float(ctaStra_NSR+ctaConv_NSR+ctaOthe_NSR+ctaNoRa_NSR)  ;;mean rainfall for all pixels
                       if ctaStra_NSR ne 0 then m_RainStrt_NSR=RTo_stra_NSR/ctaStra_NSR else m_RainStrt_NSR=0. ;;mean stratiform rain
                       if ctaConv_NSR ne 0 then m_RainConv_NSR=RTo_conv_NSR/ctaConv_NSR else m_RainConv_NSR=0. ;;mean convective rain
                       if ctaStra_NSR ne 0 or ctaConv_NSR ne 0 then $
                          m_RainCSs_NSR=total_RainCSs_NSR/float(ctaStra_NSR+ctaConv_NSR) else m_RainCSs_NSR=-999.
                       ;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       vol_Rain_All_NSR=total_RainAll_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Str_NSR=RTo_stra_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Con_NSR=RTo_conv_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       tmp_statsRain_NSR=[m_RainAll_NSR,m_RainStrt_NSR,m_RainConv_NSR,vol_Rain_All_NSR,vol_Rain_Str_NSR,vol_Rain_Con_NSR]

                       ;;*********************************************************************************************
                       ;;Now only for the Convective Core
                       ;;*******************************************************************************************
                       ;; ADDED _CV TO THESE VAR NAMES FOR CLARITY  ;; SRB 8/3/2018
                       ;; Z-R Method
                       ;;RTo_conv_R11_CV=0. & ctaConv_R11_CV=0l ;;total convective rain                                   ;; SRB
                       ;;RTo_stra_R11_CV=0. & ctaStra_R11_CV=0l ;;total stratiform rain                                   ;; SRB
                       ;;RTo_othe_R11_CV=0. & ctaOthe_R11_CV=0l ;;total other rain                                        ;; SRB
                       ;;RTo_noRa_R11_CV=0. & ctaNoRa_R11_CV=0l ;;total no_rain rain                                      ;; SRB

                       ;; ADDED _CV TO THESE VAR NAMES FOR CLARITY  ;; SRB 8/3/2018
                       ;; Near Surface Rain
                       RTo_conv_NSR_CV=0. & ctaConv_NSR_CV=0l ;;total convective rain
                       RTo_stra_NSR_CV=0. & ctaStra_NSR_CV=0l ;;total stratiform rain
                       RTo_othe_NSR_CV=0. & ctaOthe_NSR_CV=0l ;;total other rain
                       RTo_noRa_NSR_CV=0. & ctaNoRa_NSR_CV=0l ;;total no_rain rain

                       ;;Count pixels in matrix matching with storm location and locate indices (CORE STORM)
                       for i=0l,pixelsumCV_CV-1l do begin   
                          col=dondeCV_CV[i] mod nlonsFull               ;;column ID of the pixel
                          fil=long(fix(dondeCV_CV[i]/float(nlonsFull))) ;;row    ID of the pixel
                   
                          nearSrfR_Org=SrfRain_FULL[col,fil]
                          
                          t_col=where(lonsC le lonsFull_sub[col],ctaC) ;;count the pixels
                          t_row=where(latsC le latsFull_sub[fil],ctaR) 
                          if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                             colCV=(reverse(t_col))[0] 
                             rowCV=(reverse(t_row))[0]  
                             if area_CV ge thr_aCV and dim_topCV ge thr_hCV then $
                                freq_Core[colCV,rowCV,2]=freq_Core[colCV,rowCV,2]+1l else $
                                   if area_CV ge thr_aCV then freq_Core[colCV,rowCV,1]=freq_Core[colCV,rowCV,1]+1l $
                                   else if dim_topCV ge thr_hCV then freq_Core[colCV,rowCV,0]=freq_Core[colCV,rowCV,0]+1l 
                             
                             ;;locate indices (x,y,z location of individual pixels within the Convective storm)
                             pila=reform(grid_storm_FULL[col,fil,*]) ;;Matching the location within the full storm containing the core
                             w_CV=where(pila ne 0,ctaHgt)
                             if ctaHgt ge 1 then begin 
                                ;;distance between lowest pixel and ground
                                id_top1=(where(topo_lon ge lonsFull_sub[col]))[0]
                                id_top2=(where(topo_lat ge latsFull_sub[fil]))[0]
 
                                ;;%we sort out all the pixels that are higher than 2.5km above the ground
                                ;;if (hgts[w_CV[0]]-float(DEM[id_top1,id_top2])/1000.) le 2.5 then begin                    ;; SRB

                                   ;;locate the pixel in the nearest fine grid cell
                                   tmp_col=(where(float(lonsF) eq lonsFull_sub[col],ctaC))[0]
                                   tmp_row=(where(float(latsF) eq latsFull_sub[fil],ctaR))[0]
                            
                                   if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                   ;;if ctaC eq 0 or ctaR eq 0 then stop
                            
                                      ;;if nearSrfR_Org ne -9999.00 and refl_3D_FULL[col,fil,w_CV[0]] ne -999.0 then begin     ;; SRB
                                      if nearSrfR_Org ne SrfRain_fillValue then begin
                                         ;;reflectivZ=10^(refl_3D_FULL[col,fil,w_CV[0]]*0.1)  ;;%convert from dBZ to Z         ;; SRB
                               
                                         ;;here I create accumulated matrices with rain rate for each type of convective system
                                         if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective
                                            rain_NSRCore[tmp_col,tmp_row,2]=rain_NSRCore[tmp_col,tmp_row,2]+nearSrfR_Org
                                            nRai_NSRCore[tmp_col,tmp_row,2]=nRai_NSRCore[tmp_col,tmp_row,2]+1l
                                            
                                         endif else begin
                                            if area_CV ge thr_aCV then begin
                                               rain_NSRCore[tmp_col,tmp_row,1]=rain_NSRCore[tmp_col,tmp_row,1]+nearSrfR_Org
                                               nRai_NSRCore[tmp_col,tmp_row,1]=nRai_NSRCore[tmp_col,tmp_row,1]+1l
                                               
                                            endif else begin
                                               if dim_topCV ge thr_hCV then begin
                                                  rain_NSRCore[tmp_col,tmp_row,0]=rain_NSRCore[tmp_col,tmp_row,0]+nearSrfR_Org
                                                  nRai_NSRCore[tmp_col,tmp_row,0]=nRai_NSRCore[tmp_col,tmp_row,0]+1l
                                                  
                                               endif 
                                            endelse
                                         endelse
                               
                                         ;;here I create accumulated rain rate vectors
                                         if raintypeFULL[col,fil] eq CONV then begin ;;convective rain
                                            RTo_conv_NSR_CV=RTo_conv_NSR_CV+nearSrfR_Org                     &  ctaConv_NSR_CV=ctaConv_NSR_CV+1l
                                            ;;RTo_conv_R11_CV=RTo_conv_R11_CV+(reflectivZ/aRRc)^(1/bRRc)       &  ctaConv_R11_CV=ctaConv_R11_CV+1l   ;; SRB
                                            
                                            ;;tmp_rain=(reflectivZ/aRRc)^(1/bRRc)                                                        ;; SRB
                                            ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective       ;; SRB
                                            ;;   rain_R11Core[tmp_col,tmp_row,2]=rain_R11Core[tmp_col,tmp_row,2]+tmp_rain                ;; SRB
                                            ;;   nRai_R11Core[tmp_col,tmp_row,2]=nRai_R11Core[tmp_col,tmp_row,2]+1l                      ;; SRB
                                            ;;endif else begin                                                                           ;; SRB
                                            ;;   if area_CV ge thr_aCV then begin                                                        ;; SRB
                                            ;;      rain_R11Core[tmp_col,tmp_row,1]=rain_R11Core[tmp_col,tmp_row,1]+tmp_rain             ;; SRB
                                            ;;      nRai_R11Core[tmp_col,tmp_row,1]=nRai_R11Core[tmp_col,tmp_row,1]+1l                   ;; SRB
                                            ;;   endif else begin                                                                        ;; SRB
                                            ;;      if dim_topCV ge thr_hCV then begin                                                   ;; SRB
                                            ;;         rain_R11Core[tmp_col,tmp_row,0]=rain_R11Core[tmp_col,tmp_row,0]+tmp_rain          ;; SRB
                                            ;;         nRai_R11Core[tmp_col,tmp_row,0]=nRai_R11Core[tmp_col,tmp_row,0]+1l                ;; SRB
                                            ;;      endif                                                                                ;; SRB
                                            ;;   endelse                                                                                 ;; SRB
                                            ;;endelse                                                                                    ;; SRB

                                         endif else begin
                                            if raintypeFULL[col,fil] eq STRA then begin ;;stratiform rain
                                               RTo_stra_NSR_CV=RTo_stra_NSR_CV+nearSrfR_Org                  &  ctaStra_NSR_CV=ctaStra_NSR_CV+1l
                                               ;;RTo_stra_R11_CV=RTo_stra_R11_CV+(reflectivZ/aRRs)^(1/bRRs)    &  ctaStra_R11_CV=ctaStra_R11_CV+1l   ;; SRB
                                               
                                               ;;tmp_rain=(reflectivZ/aRRs)^(1/bRRs)                                                     ;; SRB
                                               ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective    ;; SRB
                                               ;;   rain_R11Core[tmp_col,tmp_row,2]=rain_R11Core[tmp_col,tmp_row,2]+tmp_rain             ;; SRB
                                               ;;   nRai_R11Core[tmp_col,tmp_row,2]=nRai_R11Core[tmp_col,tmp_row,2]+1l                   ;; SRB
                                               ;;endif else begin                                                                        ;; SRB
                                               ;;   if area_CV ge thr_aCV then begin                                                     ;; SRB
                                               ;;      rain_R11Core[tmp_col,tmp_row,1]=rain_R11Core[tmp_col,tmp_row,1]+tmp_rain          ;; SRB
                                               ;;      nRai_R11Core[tmp_col,tmp_row,1]=nRai_R11Core[tmp_col,tmp_row,1]+1l                ;; SRB
                                               ;;   endif else begin                                                                     ;; SRB
                                               ;;      if dim_topCV ge thr_hCV then begin                                                ;; SRB
                                               ;;         rain_R11Core[tmp_col,tmp_row,0]=rain_R11Core[tmp_col,tmp_row,0]+tmp_rain       ;; SRB
                                               ;;         nRai_R11Core[tmp_col,tmp_row,0]=nRai_R11Core[tmp_col,tmp_row,0]+1l             ;; SRB
                                               ;;      endif                                                                             ;; SRB
                                               ;;   endelse                                                                              ;; SRB
                                               ;;endelse                                                                                 ;; SRB
                                     
                                            endif else begin
                                               if raintypeFULL[col,fil] ge OTHER then begin ;;other rain
                                                  RTo_othe_NSR_CV=RTo_othe_NSR_CV+nearSrfR_Org               &  ctaOthe_NSR_CV=ctaOthe_NSR_CV+1l
                                                  ;;RTo_othe_R11_CV=RTo_othe_R11_CV+(reflectivZ/aRR)^(1/bRR)   &  ctaOthe_R11_CV=ctaOthe_R11_CV+1l   ;; SRB
                                                  
                                                  ;;tmp_rain=(reflectivZ/aRR)^(1/bRR)                                                    ;; SRB
                                                  ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective ;; SRB
                                                  ;;   rain_R11Core[tmp_col,tmp_row,2]=rain_R11Core[tmp_col,tmp_row,2]+tmp_rain          ;; SRB
                                                  ;;   nRai_R11Core[tmp_col,tmp_row,2]=nRai_R11Core[tmp_col,tmp_row,2]+1l                ;; SRB
                                                  ;;endif else begin                                                                     ;; SRB
                                                  ;;   if area_CV ge thr_aCV then begin                                                  ;; SRB
                                                  ;;      rain_R11Core[tmp_col,tmp_row,1]=rain_R11Core[tmp_col,tmp_row,1]+tmp_rain       ;; SRB
                                                  ;;      nRai_R11Core[tmp_col,tmp_row,1]=nRai_R11Core[tmp_col,tmp_row,1]+1l             ;; SRB
                                                  ;;   endif else begin                                                                  ;; SRB
                                                  ;;      if dim_topCV ge thr_hCV then begin                                             ;; SRB
                                                  ;;         rain_R11Core[tmp_col,tmp_row,0]=rain_R11Core[tmp_col,tmp_row,0]+tmp_rain    ;; SRB
                                                  ;;         nRai_R11Core[tmp_col,tmp_row,0]=nRai_R11Core[tmp_col,tmp_row,0]+1l          ;; SRB
                                                  ;;      endif                                                                          ;; SRB
                                                  ;;   endelse                                                                           ;; SRB
                                                  ;;endelse                                                                              ;; SRB
                                        
                                               endif else begin
                                                  if raintype[col,fil] eq raintype_noRainValue then begin ;;No rain
                                                     RTo_noRa_NSR_CV=RTo_noRa_NSR_CV+0.   &   ctaNoRa_NSR_CV=ctaNoRa_NSR_CV+1l
                                                     ;;RTo_noRa_R11_CV=RTo_noRa_R11_CV+0.   &   ctaNoRa_R11_CV=ctaNoRa_R11_CV+1l            ;; SRB
                                                     
                                                     ;;if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;deep and wide convective ;; SRB
                                                     ;;   rain_R11Core[tmp_col,tmp_row,2]=rain_R11Core[tmp_col,tmp_row,2]+0                 ;; SRB
                                                     ;;   nRai_R11Core[tmp_col,tmp_row,2]=nRai_R11Core[tmp_col,tmp_row,2]+1l                ;; SRB
                                                     ;;endif else begin                                                                     ;; SRB
                                                     ;;   if area_CV ge thr_aCV then begin                                                  ;; SRB
                                                     ;;      rain_R11Core[tmp_col,tmp_row,1]=rain_R11Core[tmp_col,tmp_row,1]+0              ;; SRB
                                                     ;;      nRai_R11Core[tmp_col,tmp_row,1]=nRai_R11Core[tmp_col,tmp_row,1]+1l             ;; SRB
                                                     ;;   endif else begin                                                                  ;; SRB
                                                     ;;      if dim_topCV ge thr_hCV then begin                                             ;; SRB
                                                     ;;         rain_R11Core[tmp_col,tmp_row,0]=rain_R11Core[tmp_col,tmp_row,0]+0           ;; SRB
                                                     ;;         nRai_R11Core[tmp_col,tmp_row,0]=nRai_R11Core[tmp_col,tmp_row,0]+1l          ;; SRB
                                                     ;;      endif                                                                          ;; SRB
                                                     ;;   endelse                                                                           ;; SRB
                                                     ;;endelse                                                                              ;; SRB
                                                  endif
                                               endelse
                                            endelse
                                         endelse
                                      endif    ;; end for missing values....
                                   endif       ;; if ctaC ne 0 and ctaR ne 0
                                ;;endif       ;; end for pixels with height > 2.5km                                                      ;; SRB
                             endif          ;;checking if there is data in the vertical
                          endif             ;;checking if there is data inside the boundaries of full storm  
                       endfor               ;;for going thry different pixelsum (indivdual pixels within storm)

                       ;;************************************************************************************************************    ;; SRB
                       ;;;;Romatschke and Houze 2011
                       ;;total_RainAll_R11=RTo_stra_R11_CV+RTo_conv_R11_CV+RTo_othe_R11_CV+RTo_noRa_R11_CV
                       ;;total_RainCSs_R11=RTo_stra_R11_CV+RTo_conv_R11_CV
                       ;;m_RainAll_R11=total_RainAll_R11/float(ctaStra_R11_CV+ctaConv_R11_CV+ctaOthe_R11_CV+ctaNoRa_R11_CV)  ;;mean rainfall for all pixels
                       ;;if ctaStra_R11_CV ne 0 then m_RainStrt_R11=RTo_stra_R11_CV/ctaStra_R11_CV else m_RainStrt_R11=0. ;;mean stratiform rain
                       ;;if ctaConv_R11_CV ne 0 then m_RainConv_R11=RTo_conv_R11_CV/ctaConv_R11_CV else m_RainConv_R11=0. ;;mean convective rain
                       ;;if ctaStra_R11_CV ne 0 or ctaConv_R11_CV ne 0 then $
                       ;;   m_RainCSs_R11=total_RainCSs_R11/float(ctaStra_R11_CV+ctaConv_R11_CV) else m_RainCSs_R11=-999.
                       ;;;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       ;;size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;vol_Rain_All_R11=total_RainAll_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Str_R11=RTo_stra_R11_CV*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Con_R11=RTo_conv_R11_CV*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       ;;tmp_statsRain_R11_CV=[m_RainAll_R11,m_RainStrt_R11,m_RainConv_R11,vol_Rain_All_R11,vol_Rain_Str_R11,vol_Rain_Con_R11]
                       
                       ;;************************************************************************************************************
                       ;;Near Surface Rain
                       total_RainAll_NSR=RTo_stra_NSR_CV+RTo_conv_NSR_CV+RTo_othe_NSR_CV+RTo_noRa_NSR_CV
                       total_RainCSs_NSR=RTo_stra_NSR_CV+RTo_conv_NSR_CV
                       m_RainAll_NSR=total_RainAll_NSR/float(ctaStra_NSR_CV+ctaConv_NSR_CV+ctaOthe_NSR_CV+ctaNoRa_NSR_CV)  ;;mean rainfall for all pixels
                       if ctaStra_NSR_CV ne 0 then m_RainStrt_NSR=RTo_stra_NSR_CV/ctaStra_NSR_CV else m_RainStrt_NSR=0. ;;mean stratiform rain
                       if ctaConv_NSR_CV ne 0 then m_RainConv_NSR=RTo_conv_NSR_CV/ctaConv_NSR_CV else m_RainConv_NSR=0. ;;mean convective rain
                       if ctaStra_NSR_CV ne 0 or ctaConv_NSR_CV ne 0 then $
                          m_RainCSs_NSR=total_RainCSs_NSR/float(ctaStra_NSR_CV+ctaConv_NSR_CV) else m_RainCSs_NSR=-999.
                       ;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       vol_Rain_All_NSR=total_RainAll_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Str_NSR=RTo_stra_NSR_CV*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Con_NSR=RTo_conv_NSR_CV*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       tmp_statsRain_NSR_CV=[m_RainAll_NSR,m_RainStrt_NSR,m_RainConv_NSR,vol_Rain_All_NSR,vol_Rain_Str_NSR,vol_Rain_Con_NSR]
                       ;;************************************************************************************************************

                       ;;*********************************************************************************************
                       ;;Here I calculate the CFAD count !!! for Full Storm!!! (only for one single Full storm)
                       if area ne areaIsFULL and cen_lon ne lonCIsFULL and cen_lat ne latCIsFull then begin
                          refl_SingleStorm=fltarr(nlonsFull,nlatsFull,nlevels)
                          refl_SingleStorm[*,*,*]=refl_3D_fillValue
                          refl_SingleStorm[w_idF]=refl_3D_FULL[w_idF]
                          if cta_Full ne npix_str[donde_Convec[ss]] then stop  ;; just to check! because this is in 3D!
                   
                          ;;here count reflectivity for each pixel that compose the storm into a matrix of CFAD
                          for i=0l,cta_Full-1l do begin  
                             col_R=where(refl_CFAD eq round(refl_3D_FULL[w_idF[i]]),ctaZ)   
                             row_H=where(alts_CFAD eq hgts_3D_FULL[w_idF[i]],ctaH)   ;;here locate the height of the pixel
                             if ctaH ne 0 and ctaZ ne 0 then $
                                if area_CV ge thr_aCV and dim_topCV ge thr_hCV then $
                                   CFAD_Full[col_R,row_H,2]=CFAD_Full[col_R,row_H,2]+1l else $
                                      if area_CV ge thr_aCV then CFAD_Full[col_R,row_H,1]=CFAD_Full[col_R,row_H,1]+1l else $
                                         if dim_topCV ge thr_hCV then CFAD_Full[col_R,row_H,0]=CFAD_Full[col_R,row_H,0]+1l  
                          endfor
                       endif

                       ;;*********************************************************************************************
                       ;;Here I calculate the CFAD count !!! for Convective component!!!
                       refl_SingleStorm=fltarr(nlonsFull,nlatsFull,nlevels)
                       refl_SingleStorm[*,*,*]=refl_3D_fillValue
                       refl_SingleStorm[w_idCV]=refl_3D_FULL[w_idCV]
                       if cta_CV ne npix_CV[donde_Convec2[ssCV]] then stop  ;; just to check! because this is in 3D!
                
                       ;;here count reflectivity for each pixel that compose the storm into a matrix of CFAD
                       for i=0l,cta_CV-1 do begin  
                          col_R=where(refl_CFAD eq round(refl_3D_FULL[w_idCV[i]]),ctaZ)   
                          row_H=where(alts_CFAD eq hgts_3D_FULL[w_idCV[i]],ctaH)   ;;here locate the height of the pixel
                          if ctaH ne 0 and ctaZ ne 0 then $
                             if area_CV ge thr_aCV and dim_topCV ge thr_hCV then $
                                CFAD_Core[col_R,row_H,2]=CFAD_Core[col_R,row_H,2]+1l else $
                                   if area_CV ge thr_aCV then CFAD_Core[col_R,row_H,1]=CFAD_Core[col_R,row_H,1]+1l else $
                                      if dim_topCV ge thr_hCV then CFAD_Core[col_R,row_H,0]=CFAD_Core[col_R,row_H,0]+1l  
                       endfor

                       areaIsFULL=area & lonCIsFULL=cen_lon & latCIsFull=cen_lat   ;;return area counter to avoid double count of fullStorm

                       if area_CV ge thr_aCV and dim_topCV ge thr_hCV then begin ;;store info ;;both Convective categories
                          ;;info_DW=[info_DW,orbit+'.'+datetime+'.'+strtrim(string(s_idCV),2)]  ;;MASK MOD
                          info_DW=[info_DW,orbit+'.'+datetime+'.'+strtrim(string(num_dwc),2)]  ;;MASK MOD
                          shape_Core_DW=[[shape_Core_DW],[tmp_shapeCV]]
                          shape_Full_DW=[[shape_Full_DW],[tmp_shape]]
                          rain_Core_DW=[[rain_Core_DW],[rain_momentCV]]
                          rain_Full_DW=[[rain_Full_DW],[rain_moment]]
                          rainTypeCore_DW=[[rainTypeCore_DW],[statsRain_CV]]
                          rainTypeFull_DW=[[rainTypeFull_DW],[statsRain]]

                          rainCore_DW_NSR=[[rainCore_DW_NSR],[tmp_statsRain_NSR_CV]]
                          rainFull_DW_NSR=[[rainFull_DW_NSR],[tmp_statsRain_NSR]]
                          ;;rainCore_DW_R11=[[rainCore_DW_R11],[tmp_statsRain_R11_CV]]                                     ;; SRB
                          ;;rainFull_DW_R11=[[rainFull_DW_R11],[tmp_statsRain_R11]]                                        ;; SRB

                       endif else begin
                          if area_CV ge thr_aCV then begin   ;;store info Wide Convective Subset (Convective area2D >= thr_aCV)
                             ;;info_WC=[info_WC,orbit+'.'+datetime+'.'+strtrim(string(s_idCV),2)]  ;;MASK MOD
                             info_WC=[info_WC,orbit+'.'+datetime+'.'+strtrim(string(num_wcc),2)]  ;;MASK MOD
                             shape_Core_WC=[[shape_Core_WC],[tmp_shapeCV]]
                             shape_Full_WC=[[shape_Full_WC],[tmp_shape]]
                             rain_Core_WC=[[rain_Core_WC],[rain_momentCV]]
                             rain_Full_WC=[[rain_Full_WC],[rain_moment]]
                             rainTypeCore_WC=[[rainTypeCore_WC],[statsRain_CV]]
                             rainTypeFull_WC=[[rainTypeFull_WC],[statsRain]]
                             
                             rainCore_WC_NSR=[[rainCore_WC_NSR],[tmp_statsRain_NSR_CV]]
                             rainFull_WC_NSR=[[rainFull_WC_NSR],[tmp_statsRain_NSR]]
                             ;;rainCore_WC_R11=[[rainCore_WC_R11],[tmp_statsRain_R11_CV]]                                  ;; SRB
                             ;;rainFull_WC_R11=[[rainFull_WC_R11],[tmp_statsRain_R11]]                                     ;; SRB

                          endif else begin
                             if dim_topCV ge thr_hCV then begin ;;store info Deep Convective Subset (Convective hgt >= thr_hCV)
                                
                                ;;info_DC=[info_DC,orbit+'.'+datetime+'.'+strtrim(string(s_idCV),2)]  ;;MASK MOD
                                info_DC=[info_DC,orbit+'.'+datetime+'.'+strtrim(string(num_dcc),2)]  ;;MASK MOD
                                shape_Core_DC=[[shape_Core_DC],[tmp_shapeCV]]
                                shape_Full_DC=[[shape_Full_DC],[tmp_shape]]
                                rain_Core_DC=[[rain_Core_DC],[rain_momentCV]]
                                rain_Full_DC=[[rain_Full_DC],[rain_moment]]
                                rainTypeCore_DC=[[rainTypeCore_DC],[statsRain_CV]]
                                rainTypeFull_DC=[[rainTypeFull_DC],[statsRain]]

                                rainCore_DC_NSR=[[rainCore_DC_NSR],[tmp_statsRain_NSR_CV]]
                                rainFull_DC_NSR=[[rainFull_DC_NSR],[tmp_statsRain_NSR]]
                                ;;rainCore_DC_R11=[[rainCore_DC_R11],[tmp_statsRain_R11_CV]]                               ;; SRB
                                ;;rainFull_DC_R11=[[rainFull_DC_R11],[tmp_statsRain_R11]]                                  ;; SRB

                             endif
                          endelse
                       endelse
                       undefine,dim_lonCV
                       undefine,dim_latCV
                       undefine,hgt_sumCV
                       undefine,grid_sum
                       undefine,donde
                       undefine,pixelsum
                       undefine,lon_sum
                       undefine,lat_sum
                       undefine,hgt_sum
                       undefine,tmp_shapeCV
                       undefine,tmp_shape
                       undefine,SrfRain_FULL
                       undefine,raintypeFULL
                       undefine,refl_3D_FULL
                       undefine,hgts_3D_FULL
                       undefine,grid_storm_FULL

                       undefine,rainCV
                       undefine,rain
                       undefine,rain_nomiss
                       undefine,stratconv
                       undefine,strats
                       undefine,convec
                       undefine,others
                       undefine,noRain
                       undefine,missin
                       undefine,pila

                       undefine,col_R
                       undefine,row_H
                       undefine,refl_SingleStorm
                       undefine,colCV
                       undefine,rowCV
                       undefine,t_col
                       undefine,t_row

                    endif      ;;endif found a storm cluster with contiguous convective pixels within theresholds 
                    undefine,s_idCV
                    undefine,w_idCV
                    undefine,singlestormgrid_CV
                    undefine,grid_sum_CV_CV
                    undefine,dondeCV_CV
                    undefine,lonCV
                    undefine,latCV
                    undefine,hgt_sumCV
                 endfor      ;;endfor loop through analyzed storms clusters that maybe are deep-wide convective
                 undefine,id_CV
                 undefine,npix_CV
                 undefine,grid_CV
                 undefine,searchNaN_CV
                 undefine,donde_Convec2
                 undefine,grid_CV
              endif      ;;endif for convective areas that could be matching the theresholds
              undefine,lonCV
              undefine,latCV
              undefine,hgt_sumCV
              undefine,dondeCV
              undefine,grid_sumCV
              undefine,singlestormgridConvective
              undefine,lonsFull_sub
              undefine,latsFull_sub
              undefine,singlestormgrid_Full
              undefine,d_latsFull
              undefine,d_lonsFull
              undefine,total_latFull
              undefine,total_lonFull
              undefine,s_idF
              undefine,w_idF
           endfor            ;;endfor loop thru convective volumes greater than 25 pixels
        endif                ;;end if the orbit does have convective and reflec > 0 ==>raining
        
        ;;---------------------
        ;;---------------------
        print,'Id-Shallow Isol'
        ;;---------------------
        ;;---------------------
        ;;here it is going to identify shallow convective part
        if cta_ShIs ne 0 then begin ;;here analyze only when shallow pixels exist
           ;;here I only choose shallow convective volumes with more than 2 pixels (same as for DCC-WCC)
           donde_shallow=where(npix_str ge 2l,ctaShallow)
       
           for ss=0l,ctaShallow-1 do begin 
              s_idF=id_storm[donde_shallow[ss]]
              w_idF=where(grid_storm eq s_idF,cta1)
              if cta1 ne npix_str[donde_shallow[ss]] then stop  ;; just to check!
              singlestormgrid=lonarr(nlons,nlats,nlevels)
              singlestormgrid[w_idF]=s_idF

              ;;Here I will subset the singleFULL storm found to reduce the size of the matrix
              total_lonFull=total(total(singlestormgrid,2),2) 
              d_lonsFull=where(total_lonFull gt 0l,nlonsFull) ;;this are the horizontal positions within a selcted contiguous area
              total_latFull=total(total(singlestormgrid,1),2) 
              d_latsFull=where(total_latFull gt 0l,nlatsFull)
         
              singlestormgrid_Full=lonarr(nlonsFull,nlatsFull,nlevels)
              singlestormgrid_Full=singlestormgrid[d_lonsFull,d_latsFull,*]
              w_idF=where(singlestormgrid_Full eq s_idF,cta_Full)
              lonsFull_sub=lons[d_lonsFull]
              latsFull_sub=lats[d_latsFull]
         
              undefine,singlestormgrid

              ;;Identify the Shallow pixels in storm (refl >= 0 and All shallow convectve pixels)
              singlestormgridShallow=singlestormgrid_Full
              singlestormgridShallow[where(refl_3D[d_lonsFull,d_latsFull,*]  lt 0.)]=0l
              singlestormgridShallow[where(rain_type3D[d_lonsFull,d_latsFull,*] ne SHIS)]=0l   ;;everything but shallow ==0l
               
              ;;2D horizontal array of number of Shallow pixels within storm
              grid_sumSH=total(singlestormgridShallow,3) 
              dondeSH=where(grid_sumSH gt 0,pixelsumSH) ;;pixelsum=number of pixels in a 2D proj
         
              if pixelsumSH ge 2 then begin  ;;shallow conv volumes with more than 2 pixels (to avoid smaller areas)
                 lonSH=total(total(singlestormgridShallow,2),2) 
                 lon_SH=(max(lonsFull_sub[where(lonSH gt 0l)])+min(lonsFull_sub[where(lonSH gt 0l)]))/2. ;;center of Shallow

                 latSH=total(total(singlestormgridShallow,1),2) 
                 lat_SH=(max(latsFull_sub[where(latSH gt 0l)])+min(latsFull_sub[where(latSH gt 0l)]))/2. ;;center of Shallow

                 size_pixels=deg2km(pixDeg,lon_SH,lat_SH)
                 area_SH=pixelsumSH*(size_pixels[0]*size_pixels[1]) ;;stratiform area in km  
              endif else begin
                 area_SH=0.                                 ;;to avoid regions with two few pixels togheter
                 size_pixels=[pixKm,pixKm]                  ;; this is set so next if loop not entered when area=0
              endelse
                 
              ;;This is original way; new way more restrictive (pixArea*2 =60.5)        ;; SRB
              ;;if area_Sh gt 40. then begin ;;  ~30km2 is the size of a single pixel (Avoid single pixels!)
              if area_Sh ge (2*(size_pixels[0]*size_pixels[1])) then begin ;;  Avoid single pixels
              ;;if area_Sh ge (pixArea*2) then begin ;;  pixArea (~30km2) is the size of a single pixel (Avoid single pixels!)
                 ;;Does the shallow isolated area belongs to a single storm? (acomplish contiguous pixel condition)
                 ;;Again, uses the shallow pixels to identify contiguous pixels
                 findStormNew,refl_3D=singlestormgridShallow,$
                              ;;refl_3D_fillValue=refl_3D_fillValue,$
                              id_storm=id_SH,$
                              npix_str=npix_SH,$
                              grid_storm=grid_SH
                 searchNaN_SH=where(grid_SH lt 0, nanCnt)
                 if nanCnt gt 0 then grid_SH[searchNaN_SH]=0l ;;%set NaNs to zero in the  grid matrix 
        
                 ;;identify only volumes with more than 4 pixels - THIS DOESN'T MAKE SENSE  ;;SRB
                 donde_Shallow2=where(npix_SH ge 2l,ctaShallow2) ;;permiting to have 2 contiguous pixels of greater

                 areaIsFULL=0. & lonCIsFULL=0. & latCIsFull=0. ;;this is to locate only ONE full storm

                 for ssSH=0,ctaShallow2-1 do begin ;;only goes thru volumes greater than 2 pixels (speed the process)
                    s_id_SH=id_SH[donde_Shallow2[ssSH]]
                    w_id_SH=where(grid_SH eq s_id_SH,cta1_SH)
                    if cta1_SH ne npix_SH[donde_Shallow2[ssSH]] then stop  ;; just to check!
              
                    singlestormgrid_SH=lonarr(nlonsFull,nlatsFull,nlevels)
                    singlestormgrid_SH[w_id_SH]=s_id_SH
              
                    ;;2D horizontal array of number of Shallow pixels within storm
                    grid_sum_SH_SH=total(singlestormgrid_SH,3) 
                    dondeSH_SH=where(grid_sum_SH_SH gt 0,pixelsumSH_SH) ;;pixelsum=number of pixels in a 2D proj
               
                    ;; this is going to be done because the calculation of rain and LH for ALL-Events is done considering
                    ;; 2 horizontal pixels or more!!
                    if pixelsumSH_SH ge 2 then begin
                       lonSH=total(total(singlestormgrid_SH,2),2) 
                       lonC_SH=(max(lonsFull_sub[where(lonSH gt 0l)])+min(lonsFull_sub[where(lonSH gt 0l)]))/2. ;;center of Shallow
                 
                       latSH=total(total(singlestormgrid_SH,1),2) 
                       latC_SH=(max(latsFull_sub[where(latSH gt 0l)])+min(latsFull_sub[where(latSH gt 0l)]))/2. ;;center of Shallow
                 
                       size_pixels=deg2km(pixDeg,lonC_SH,latC_SH)
                       area_SH=pixelsumSH_SH*(size_pixels[0]*size_pixels[1]) ;;Broad Shallow area in km  
                    endif else begin
                       area_SH=0.
                       size_pixels=[pixKm,pixKm]  ;; this is set so next if loop not entered when area=0
                    endelse 

                    ;;now it tries to identify the real Shallow contiguous area
                    ;;if area_SH ge 30. then begin ;; This is original way        ;; SRB
                    if area_SH ge (2*(size_pixels[0]*size_pixels[1])) then begin ;; area of one pixel is pixArea (~30km^2)
                       ;;calculate the dimensions for Shallow pixels in the cluster
                       ;;******************************************************************************************
                       dim_lonSH=(max(lonsFull_sub[where(lonSH gt 0l)])-min(lonsFull_sub[where(lonSH gt 0l)]))+pixDeg
                       dim_latSH=(max(latsFull_sub[where(latSH gt 0l)])-min(latsFull_sub[where(latSH gt 0l)]))+pixDeg

                       hgt_sumSH=total(total(singlestormgrid_SH,2),1)
                       dim_hgtSH=max(hgts[where(hgt_sumSH gt 0l)])-min(hgts[where(hgt_sumSH gt 0l)])
                       dim_topSH=max(hgts[where(hgt_sumSH gt 0l)])                                   ;;###############
                       dim_botSH=min(hgts[where(hgt_sumSH gt 0l)])                                   ;;###############
                
                       ;;calculates the elevation of the terrain for the center of the storm  !21
                       id_top1=(where(topo_lon ge lonC_SH))[0] & id_top2=(where(topo_lat ge latC_SH))[0]
                       terr_hgtSH=DEM[id_top1,id_top2]                              ;;elevation in meters
                       if terr_hgtSH eq 0 then land_oceanSH=0 else land_oceanSH=1   ;;ocean=0 or land=1!

                       tmp_shapeSH=[lonC_SH,latC_SH,area_SH,dim_topSH,dim_botSH,dim_lonSH,dim_latSH,terr_hgtSH,land_oceanSH]

                       ;;num_shi = num_shi++
                       ;; add storm id to shi_mask if makeCoreMasks is set
                       ;;if makeCoreMasks then begin
                          ;;shi_mask_sub = shi_mask[d_lonsFull,d_latsFull,0]
                          ;;shi_mask_sub[dondeSH_SH]=num_shi
                          ;;shi_mask[d_lonsFull,d_latsFull,0] = shi_mask_sub
                          ;;undefine,shi_mask_sub
                       ;;endif

                       ;;now it calculate the dimensions for full storm (all pixels in the cluster)
                       ;;******************************************************************************************
                       grid_sum=total(singlestormgrid_Full,3) ;;2D horiz array of number of all pixels within storm
                       donde=where(grid_sum gt 0,pixelsum)    ;;pixelsum=number of pixels in a 2D proj
                
                       lon_sum=total(total(singlestormgrid_Full,2),2)   
                       cen_lon=(max(lonsFull_sub[where(lon_sum gt 0l)])+min(lonsFull_sub[where(lon_sum gt 0l)]))/2.     ;;longitude coordinate
                       dim_lon=(max(lonsFull_sub[where(lon_sum gt 0l)])-min(lonsFull_sub[where(lon_sum gt 0l)]))+pixDeg ;;longitudinal extents
                
                       lat_sum=total(total(singlestormgrid_Full,1),2)   
                       cen_lat=(max(latsFull_sub[where(lat_sum gt 0l)])+min(latsFull_sub[where(lat_sum gt 0l)]))/2.     ;;latitude coordinate
                       dim_lat=(max(latsFull_sub[where(lat_sum gt 0l)])-min(latsFull_sub[where(lat_sum gt 0l)]))+pixDeg ;;latitudinal extents
                
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       area=pixelsum*size_pixels[0]*size_pixels[1] ;;in km horizontal area of selected storm
                
                       hgt_sum=total(total(singlestormgrid_Full,2),1)
                       dim_hgt=(max(hgts[where(hgt_sum gt 0l)])-min(hgts[where(hgt_sum gt 0l)])) ;;%vertical extent in km
                       dim_top=max(hgts[where(hgt_sum gt 0l)])                                   ;;###############          
                       dim_bot=min(hgts[where(hgt_sum gt 0l)])                                   ;;###############
                
                       ;;%calculates the elevation of the terrain for the center of the storm  !21
                       id_top1=(where(topo_lon ge cen_lon))[0] & id_top2=(where(topo_lat ge cen_lat))[0]
                       terr_hgt=DEM[id_top1,id_top2]                          ;;elevation in meters
                       if terr_hgt eq 0 then land_ocean=0 else land_ocean=1   ;;mask  ocean=0 or land=1!

                       tmp_shape=[cen_lon,cen_lat,area,dim_top,dim_bot,dim_lon,dim_lat,terr_hgt,land_ocean]

                       ;;******************************************************************************************
                       ;;Statistics for rain types and rain rates *************************************************
                       ;;Statistics over Convective subset!
                       SrfRain_FULL=SrfRain[d_lonsFull,d_latsFull,*]
                       raintypeFULL=raintype[d_lonsFull,d_latsFull,*]
                       refl_3D_FULL=refl_3D[d_lonsFull,d_latsFull,*]
                       hgts_3D_FULL=hgts_3D[d_lonsFull,d_latsFull,*]

                       stratconv=raintypeFULL[dondeSH_SH]                     ;;type of rain in each 2D pixel that compose the storm
                       strats=where(stratconv eq STRA,ctaStr)                 ;;stratiform
                       convec=where(stratconv eq CONV,ctaCon)                 ;;convective
                       others=where(stratconv ge OTHER,ctaOth)                ;;other type
                       noRain=where(stratconv eq raintype_noRainValue,ctaNoR) ;;no rain
                       missin=where(stratconv eq raintype_fillValue,ctaMis)   ;;missing value

                       rainSH=SrfRain_FULL[dondeSH_SH] ;;this is based on Near Surface Rain
                       rain_nomiss=where(rainSH ne SrfRain_fillValue,Rmiss)
                       
                       ;;here I calculate moments of simple rain within storm
                       if Rmiss ge 2 then rain_momentSH=[mean(rainSH[rain_nomiss]),stdev(rainSH[rain_nomiss]),$
                                                         max(rainSH[rain_nomiss]),min(rainSH[rain_nomiss]),float(Rmiss),$
                                                         float(ctaStr),float(ctaCon)] $
                       else rain_momentSH=[-9999.,-9999.,-9999.,-9999.,-9999.,-9999.,-9999.]
                
                       ;;here I calculate rainrate sums in [mm/hr]
                       if ctaStr ne 0 then RTo_stra=total(rainSH[strats]) else RTo_stra=0.   ;;total stratiform rain
                       if ctaCon ne 0 then RTo_conv=total(rainSH[convec]) else RTo_conv=0.   ;;total convective rain
                       if ctaOth ne 0 then RTo_othe=total(rainSH[others]) else RTo_othe=0.   ;;total other rain
                       if ctaNoR ne 0 then RTo_noRa=total(rainSH[noRain]) else RTo_noRa=0.   ;;total no_rain rain
                       total_RainAll=RTo_stra+RTo_conv+RTo_othe+RTo_noRa
                       total_RainCSs=RTo_stra+RTo_conv

                       m_RainAll=total_RainAll/(ctaStr+ctaCon+ctaOth+ctaNoR+ctaMis)        ;;mean rainfall all pixels
                       if ctaStr ne 0 then m_RainStrt=RTo_stra/ctaStr else m_RainStrt=0.   ;;mean stratiform rain
                       if ctaCon ne 0 then m_RainConv=RTo_conv/ctaCon else m_RainConv=0.   ;;mean convective rain
                       if ctaStr ne 0 or ctaCon ne 0 then m_RainCSs=total_RainCSs/(ctaStr+ctaCon) $
                       else m_RainCSs=-999.
                
                       ;;here calculate volumetric values of rain  ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,lonC_SH,latC_SH)
                       vol_Rain_All=total_RainAll*(size_pixels[0]*size_pixels[1])/secsPerHr  
                       vol_Rain_Str=RTo_stra*(size_pixels[0]*size_pixels[1])/secsPerHr 
                       vol_Rain_Con=RTo_conv*(size_pixels[0]*size_pixels[1])/secsPerHr
   
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       statsRain_SH=[m_RainAll,m_RainStrt,m_RainConv,vol_Rain_All,vol_Rain_Str,vol_Rain_Con]

                       ;;Satisitics for rain types and rain rates **********************************************
                       ;;***************************************************************************************
                       stratconv=raintypeFULL[donde]                          ;;type of rain in each 2D pixel that compose the storm
                       strats=where(stratconv eq STRA,ctaStr)                 ;;stratiform
                       convec=where(stratconv eq CONV,ctaCon)                 ;;convective
                       others=where(stratconv ge OTHER,ctaOth)                ;;other type  Including shallow rain!
                       noRain=where(stratconv eq raintype_noRainValue,ctaNoR) ;;no rain
                       missin=where(stratconv eq raintype_fillValue,ctaMis)   ;;missing value
                  
                       rain=SrfRain_FULL[donde] ;;this is based on Near Surface Rain
                       rain_nomiss=where(rain ne SrfRain_fillValue,Rmiss)
                       
                       ;;here I calculate moments of simple rain within storm
                       if Rmiss ge 2 then rain_moment=[mean(rain[rain_nomiss]),stdev(rain[rain_nomiss]),$
                                                       max(rain[rain_nomiss]),min(rain[rain_nomiss]),float(Rmiss), $
                                                       float(ctaStr),float(ctaCon)] $
                       else rain_moment=[-9999.,-9999.,-9999.,-9999.,-9999.,-9999.,-9999.]
                  
                       ;;here I calculate rainrate sums in [mm/hr]
                       if ctaStr ne 0 then RTo_stra=total(rain[strats]) else RTo_stra=0.   ;;total stratiform rain
                       if ctaCon ne 0 then RTo_conv=total(rain[convec]) else RTo_conv=0.   ;;total convective rain
                       if ctaOth ne 0 then RTo_othe=total(rain[others]) else RTo_othe=0.   ;;total other rain
                       if ctaNoR ne 0 then RTo_noRa=total(rain[noRain]) else RTo_noRa=0.   ;;total no_rain rain
                       total_RainAll=RTo_stra+RTo_conv+RTo_othe+RTo_noRa
                       total_RainCSs=RTo_stra+RTo_conv
                  
                       m_RainAll=total_RainAll/(ctaStr+ctaCon+ctaOth+ctaNoR+ctaMis)        ;;mean rainfall all pixels
                       if ctaStr ne 0 then m_RainStrt=RTo_stra/ctaStr else m_RainStrt=0.   ;;mean stratiform rain
                       if ctaCon ne 0 then m_RainConv=RTo_conv/ctaCon else m_RainConv=0.   ;;mean convective rain
                       if ctaStr ne 0 or ctaCon ne 0 then m_RainCSs=total_RainCSs/(ctaStr+ctaCon) $
                       else m_RainCSs=-999.
                  
                       ;;here calculate volumetric values of rain  ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;area=pixelsum*size_pixels[0]*size_pixels[1] ;;in km horizontal area of selected storm

                       vol_Rain_All=total_RainAll*(size_pixels[0]*size_pixels[1])/secsPerHr  
                       vol_Rain_Str=RTo_stra*(size_pixels[0]*size_pixels[1])/secsPerHr 
                       vol_Rain_Con=RTo_conv*(size_pixels[0]*size_pixels[1])/secsPerHr

                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       statsRain=[m_RainAll,m_RainStrt,m_RainConv,vol_Rain_All,vol_Rain_Str,vol_Rain_Con]

                       ;;**********************************************************************************************
                       grid_storm_FULL=lonarr(nlonsFull,nlatsFull,nlevels)
                       grid_storm_FULL=grid_storm[d_lonsFull,d_latsFull,*]

                       ;;----------------------------------------------------------------------------------------------
                       ;;MOVED THIS CODE INTO IF LOOP WHERE IT NEEDS TO BE TO GET RID OF NaN OUTPUTS
                       ;;Statistics for rain rates Calculated using Z-R method 
                       ;;we need to get all the needed information for the lowest pixels in each storm
                       ;;locate the indices of the lowest pixels with data   (Full Storm) (Romatschke and Houze 2011
                       ;;;;RTo_conv_R11=0. & ctaConv_R11=0l   ;;total convective rain              ;; SRB
                       ;;;;RTo_stra_R11=0. & ctaStra_R11=0l   ;;total stratiform rain              ;; SRB
                       ;;;;RTo_othe_R11=0. & ctaOthe_R11=0l   ;;total other rain                   ;; SRB
                       ;;;;RTo_noRa_R11=0. & ctaNoRa_R11=0l   ;;total no_rain rain                 ;; SRB
                
                       ;;Near Surface Rain
                       ;;RTo_conv_NSR=0. & ctaConv_NSR=0l   ;;total convective rain
                       ;;RTo_stra_NSR=0. & ctaStra_NSR=0l   ;;total stratiform rain
                       ;;RTo_othe_NSR=0. & ctaOthe_NSR=0l   ;;total other rain
                       ;;RTo_noRa_NSR=0. & ctaNoRa_NSR=0l   ;;total no_rain rain
                       ;;----------------------------------------------------------------------------------------------

                       ;;print,'BEFORE STORM CALCS: ss=',string(strtrim(ss,2)),' and ssSH=',string(strtrim(ssSH,2)), $
                       ;;      ' and area=',string(strtrim(area,2)),' and areaIsFULL=',string(strtrim(areaIsFULL,2))

                       ;;here also is accounted n_pixels in matrix matching with storm location and locate indices (lowest pixels           
                       ;;this is only evaluated ONCE! (accouting for a single FULL storm)
                       if area ne areaIsFULL and cen_lon ne lonCIsFULL and cen_lat ne latCIsFull then begin

                          ;;print,'     RUNNING THROUGH CALCS FOR STORM - SHOULD BE ONLY ONCE'
                          
                          ;;Statistics for rain rates Calculated using Z-R method 
                          ;;we need to get all the needed information for the lowest pixels in each storm
                          ;;locate the indices of the lowest pixels with data   (Full Storm) (Romatschke and Houze 2011
                          ;;RTo_conv_R11=0. & ctaConv_R11=0l   ;;total convective rain              ;; SRB
                          ;;RTo_stra_R11=0. & ctaStra_R11=0l   ;;total stratiform rain              ;; SRB
                          ;;RTo_othe_R11=0. & ctaOthe_R11=0l   ;;total other rain                   ;; SRB
                          ;;RTo_noRa_R11=0. & ctaNoRa_R11=0l   ;;total no_rain rain                 ;; SRB
                
                          ;;Near Surface Rain
                          RTo_conv_NSR=0. & ctaConv_NSR=0l   ;;total convective rain
                          RTo_stra_NSR=0. & ctaStra_NSR=0l   ;;total stratiform rain
                          RTo_othe_NSR=0. & ctaOthe_NSR=0l   ;;total other rain
                          RTo_noRa_NSR=0. & ctaNoRa_NSR=0l   ;;total no_rain rain

                          for i=0l,pixelsum-1l do begin   
                             col=donde[i] mod nlonsFull               ;;column ID of the pixel
                             fil=long(fix(donde[i]/float(nlonsFull))) ;;row    ID of the pixel
                      
                             nearSrfR_Org=SrfRain_FULL[col,fil]
                      
                             t_col=where(lonsC le lonsFull_sub[col],ctaC) ;;count the pixels
                             t_row=where(latsC le latsFull_sub[fil],ctaR) 
                             if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                colCV=(reverse(t_col))[0] 
                                rowCV=(reverse(t_row))[0] 
                                freq_Full[colCV,rowCV,4]=freq_Full[colCV,rowCV,4]+1l 

                                ;;locate indices (x,y,z location of individual pixels within the FULL storm)
                                pila=reform(grid_storm_FULL[col,fil,*])
                                w_CV=where(pila ne 0,ctaHgt)
                                if ctaHgt ge 1 then begin 
                                   ;;distance between lowest pixel and ground
                                   id_top1=(where(topo_lon ge lonsFull_sub[col]))[0]
                                   id_top2=(where(topo_lat ge latsFull_sub[fil]))[0]

                                   ;;%we sort out all the pixels that are higher than 2.5km above the ground
                                   ;;if (hgts[w_CV[0]]-float(DEM[id_top1,id_top2])/1000.) le 2.5 then begin                  ;; SRB
                                      ;;locate the pixel in the nearest fine grid cell
                                      tmp_col=(where(float(lonsF) eq lonsFull_sub[col],ctaC))[0]
                                      tmp_row=(where(float(latsF) eq latsFull_sub[fil],ctaR))[0]
                                      
                                      if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                      ;;if ctaC eq 0 or ctaR eq 0 then stop
                               
                                         ;;if nearSrfR_Org ne -9999.00 and refl_3D_FULL[col,fil,w_CV[0]] ne -999.0 then begin   ;; SRB
                                         if nearSrfR_Org ne SrfRain_fillValue then begin
                                            ;;reflectivZ=10^(refl_3D_FULL[col,fil,w_CV[0]]*0.1)  ;;%convert from dBZ to Z
                                  
                                            ;;here I create accumulated matrices with rain rate Near Surf Rain
                                            rain_NSRFull[tmp_col,tmp_row,4]=rain_NSRFull[tmp_col,tmp_row,4]+nearSrfR_Org
                                            nRai_NSRFull[tmp_col,tmp_row,4]=nRai_NSRFull[tmp_col,tmp_row,4]+1l
                                  
                                            ;;here I create accumulated rain rate vectors with R11 method
                                            if raintypeFULL[col,fil] eq CONV then begin ;;convective rain
                                               RTo_conv_NSR=RTo_conv_NSR+nearSrfR_Org                     &  ctaConv_NSR=ctaConv_NSR+1l
                                               ;;RTo_conv_R11=RTo_conv_R11+(reflectivZ/aRRc)^(1/bRRc)       &  ctaConv_R11=ctaConv_R11+1l ;; SRB

                                               ;;tmp_rain=(reflectivZ/aRRc)^(1/bRRc)                                                      ;; SRB
                                               ;;rain_R11Full[tmp_col,tmp_row,4]=rain_R11Full[tmp_col,tmp_row,4]+tmp_rain                 ;; SRB
                                               ;;nRai_R11Full[tmp_col,tmp_row,4]=nRai_R11Full[tmp_col,tmp_row,4]+1l                       ;; SRB
                                            endif else begin
                                               if raintypeFULL[col,fil] eq STRA then begin ;;stratiform rain
                                                  RTo_stra_NSR=RTo_stra_NSR+nearSrfR_Org                  &  ctaStra_NSR=ctaStra_NSR+1l
                                                  ;;RTo_stra_R11=RTo_stra_R11+(reflectivZ/aRRs)^(1/bRRs)    &  ctaStra_R11=ctaStra_R11+1l ;; SRB

                                                  ;;tmp_rain=(reflectivZ/aRRs)^(1/bRRs)                                                   ;; SRB
                                                  ;;rain_R11Full[tmp_col,tmp_row,4]=rain_R11Full[tmp_col,tmp_row,4]+tmp_rain              ;; SRB
                                                  ;;nRai_R11Full[tmp_col,tmp_row,4]=nRai_R11Full[tmp_col,tmp_row,4]+1l                    ;; SRB
                                               endif else begin
                                                  if raintypeFULL[col,fil] ge OTHER then begin ;;other rain
                                                     RTo_othe_NSR=RTo_othe_NSR+nearSrfR_Org               &  ctaOthe_NSR=ctaOthe_NSR+1l
                                                     ;;RTo_othe_R11=RTo_othe_R11+(reflectivZ/aRR)^(1/bRR)   &  ctaOthe_R11=ctaOthe_R11+1l ;; SRB

                                                     ;;tmp_rain=(reflectivZ/aRR)^(1/bRR)                                                  ;; SRB
                                                     ;;rain_R11Full[tmp_col,tmp_row,4]=rain_R11Full[tmp_col,tmp_row,4]+tmp_rain           ;; SRB
                                                     ;;nRai_R11Full[tmp_col,tmp_row,4]=nRai_R11Full[tmp_col,tmp_row,4]+1l                 ;; SRB
                                                  endif else begin
                                                     if raintype[col,fil] eq raintype_noRainValue then begin ;;No rain
                                                        RTo_noRa_NSR=RTo_noRa_NSR+0.   &   ctaNoRa_NSR=ctaNoRa_NSR+1l
                                                        ;;RTo_noRa_R11=RTo_noRa_R11+0.   &   ctaNoRa_R11=ctaNoRa_R11+1l                   ;; SRB

                                                        ;;rain_R11Full[tmp_col,tmp_row,4]=rain_R11Full[tmp_col,tmp_row,4]+0               ;; SRB
                                                        ;;nRai_R11Full[tmp_col,tmp_row,4]=nRai_R11Full[tmp_col,tmp_row,4]+1l              ;; SRB
                                                     endif
                                                  endelse
                                               endelse
                                            endelse
                                         endif    ;; end for missing values....
                                      endif       ;; if ctaC ne 0 and ctaR ne 0
                                   ;;endif        ;; end for pixels with height > 2.5km
                                endif
                             endif
                          endfor
                       endif  ;;end for flag of accounting for a single FULL storm...
                       
                       ;;print,'AT BOTTOM OF STORM CALCS: ss=',string(strtrim(ss,2)),' and ssSH=',string(strtrim(ssSH,2)), $
                       ;;      ' and area=',string(strtrim(area,2)),' and areaIsFull=',string(strtrim(areaIsFull,2))
                       ;;;;print,'     ctaStra_R11=',string(strtrim(ctaStra_R11,2)),' and ctaConv_R11=',string(strtrim(ctaConv_R11,2)), $
                       ;;;;     ' ctaOthe_R11=',string(strtrim(ctaOthe_R11,2)),' ctaNoRa_R11=',string(strtrim(ctaNoRa_R11,2))
                       ;;print,'     ctaStra_NSR=',string(strtrim(ctaStra_NSR,2)),' and ctaConv_NSR=',string(strtrim(ctaConv_NSR,2)), $
                       ;;     ' ctaOthe_NSR=',string(strtrim(ctaOthe_NSR,2)),' ctaNoRa_NSR=',string(strtrim(ctaNoRa_NSR,2))
                       
                       ;;************************************************************************************************************  ;; SRB
                       ;;Romatschke and Houze 2011
                       ;;total_RainAll_R11=RTo_stra_R11+RTo_conv_R11+RTo_othe_R11+RTo_noRa_R11
                       ;;total_RainCSs_R11=RTo_stra_R11+RTo_conv_R11
                       ;;m_RainAll_R11=total_RainAll_R11/float(ctaStra_R11+ctaConv_R11+ctaOthe_R11+ctaNoRa_R11)  ;;mean rainfall for all pixels
                       ;;if ctaStra_R11 ne 0 then m_RainStrt_R11=RTo_stra_R11/ctaStra_R11 else m_RainStrt_R11=0. ;;mean stratiform rain
                       ;;if ctaConv_R11 ne 0 then m_RainConv_R11=RTo_conv_R11/ctaConv_R11 else m_RainConv_R11=0. ;;mean convective rain
                       ;;if ctaStra_R11 ne 0 or ctaConv_R11 ne 0 then $
                       ;;   m_RainCSs_R11=total_RainCSs_R11/float(ctaStra_R11+ctaConv_R11) else m_RainCSs_R11=-999.
                       ;;;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       ;;size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;vol_Rain_All_R11=total_RainAll_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Str_R11=RTo_stra_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Con_R11=RTo_conv_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       ;;tmp_statsRain_R11=[m_RainAll_R11,m_RainStrt_R11,m_RainConv_R11,vol_Rain_All_R11,vol_Rain_Str_R11,vol_Rain_Con_R11]

                       ;;************************************************************************************************************
                       ;;Near Surface Rain
                       total_RainAll_NSR=RTo_stra_NSR+RTo_conv_NSR+RTo_othe_NSR+RTo_noRa_NSR
                       total_RainCSs_NSR=RTo_stra_NSR+RTo_conv_NSR
                       m_RainAll_NSR=total_RainAll_NSR/float(ctaStra_NSR+ctaConv_NSR+ctaOthe_NSR+ctaNoRa_NSR)  ;;mean rainfall for all pixels
                       if ctaStra_NSR ne 0 then m_RainStrt_NSR=RTo_stra_NSR/ctaStra_NSR else m_RainStrt_NSR=0. ;;mean stratiform rain
                       if ctaConv_NSR ne 0 then m_RainConv_NSR=RTo_conv_NSR/ctaConv_NSR else m_RainConv_NSR=0. ;;mean convective rain
                       if ctaStra_NSR ne 0 or ctaConv_NSR ne 0 then $
                          m_RainCSs_NSR=total_RainCSs_NSR/float(ctaStra_NSR+ctaConv_NSR) else m_RainCSs_NSR=-999.
                       ;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       vol_Rain_All_NSR=total_RainAll_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Str_NSR=RTo_stra_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Con_NSR=RTo_conv_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       tmp_statsRain_NSR=[m_RainAll_NSR,m_RainStrt_NSR,m_RainConv_NSR,vol_Rain_All_NSR,vol_Rain_Str_NSR,vol_Rain_Con_NSR]

                       ;;*********************************************************************************************
                       ;;Now only for the Convective Core
                       ;;*******************************************************************************************
                       ;; ADDED _SH TO THESE VAR NAMES FOR CLARITY  ;; SRB 8/3/2018
                       ;; Z-R Method
                       ;;RTo_conv_R11_SH=0. & ctaConv_R11_SH=0l ;;total convective rain                                              ;; SRB
                       ;;RTo_stra_R11_SH=0. & ctaStra_R11_SH=0l ;;total stratiform rain                                              ;; SRB
                       ;;RTo_othe_R11_SH=0. & ctaOthe_R11_SH=0l ;;total other rain                                                   ;; SRB
                       ;;RTo_noRa_R11_SH=0. & ctaNoRa_R11_SH=0l ;;total no_rain rain                                                 ;; SRB

                       ;; ADDED _SH TO THESE VAR NAMES FOR CLARITY  ;; SRB 8/3/2018
                       ;; Near Surface Rain
                       RTo_conv_NSR_SH=0. & ctaConv_NSR_SH=0l ;;total convective rain
                       RTo_stra_NSR_SH=0. & ctaStra_NSR_SH=0l ;;total stratiform rain
                       RTo_othe_NSR_SH=0. & ctaOthe_NSR_SH=0l ;;total other rain
                       RTo_noRa_NSR_SH=0. & ctaNoRa_NSR_SH=0l ;;total no_rain rain

                       ;;Count pixels in matrix matching with storm location and locate indices (CORE STORM)
                       for i=0l,pixelsumSH_SH-1l do begin   
                          col=dondeSH_SH[i] mod nlonsFull               ;;column ID of the pixel
                          fil=long(fix(dondeSH_SH[i]/float(nlonsFull))) ;;row    ID of the pixel
                   
                          nearSrfR_Org=SrfRain_FULL[col,fil]
 
                          t_col=where(lonsC le lonsFull_sub[col],ctaC) ;;count the pixels
                          t_row=where(latsC le latsFull_sub[fil],ctaR) 
                          if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                             colCV=(reverse(t_col))[0] 
                             rowCV=(reverse(t_row))[0]  

                             freq_Core[colCV,rowCV,4]=freq_Core[colCV,rowCV,4]+1l

                             ;;locate indices (x,y,z location of individual pixels within the Convective storm)
                             pila=reform(grid_storm_FULL[col,fil,*]) ;;Matching the location within the full storm containing the core
                             w_CV=where(pila ne 0,ctaHgt)
                             if ctaHgt ge 1 then begin 
                                ;;distance between lowest pixel and ground
                                id_top1=(where(topo_lon ge lonsFull_sub[col]))[0]
                                id_top2=(where(topo_lat ge latsFull_sub[fil]))[0]
 
                                ;;%we sort out all the pixels that are higher than 2.5km above the ground
                                ;;if (hgts[w_CV[0]]-float(DEM[id_top1,id_top2])/1000.) le 2.5 then begin                                 ;; SRB

                                   ;;locate the pixel in the nearest fine grid cell
                                   tmp_col=(where(float(lonsF) eq lonsFull_sub[col],ctaC))[0]
                                   tmp_row=(where(float(latsF) eq latsFull_sub[fil],ctaR))[0]
                            
                                   if ctaC ne 0 and ctaR ne 0 then begin ;;just to make sure pixel located within boundaries
                                   ;;if ctaC eq 0 or ctaR eq 0 then stop
                            
                                      ;;if nearSrfR_Org ne -9999.00 and refl_3D_FULL[col,fil,w_CV[0]] ne -999.0 then begin                  ;; SRB
                                      if nearSrfR_Org ne SrfRain_fillValue then begin
                                         ;;reflectivZ=10^(refl_3D_FULL[col,fil,w_CV[0]]*0.1)  ;;%convert from dBZ to Z                      ;; SRB
                               
                                         ;;here I create accumulated matrices with rain rate for shallow Near Surf Rain
                                         rain_NSRCore[tmp_col,tmp_row,4]=rain_NSRCore[tmp_col,tmp_row,4]+nearSrfR_Org
                                         nRai_NSRCore[tmp_col,tmp_row,4]=nRai_NSRCore[tmp_col,tmp_row,4]+1l

                                         ;;here I create accumulated rain rate vectors
                                         if raintypeFULL[col,fil] eq CONV then begin ;;convective rain
                                            RTo_conv_NSR_SH=RTo_conv_NSR_SH+nearSrfR_Org                     &  ctaConv_NSR_SH=ctaConv_NSR_SH+1l
                                            ;;RTo_conv_R11_SH=RTo_conv_R11_SH+(reflectivZ/aRRc)^(1/bRRc)       &  ctaConv_R11_SH=ctaConv_R11_SH+1l      ;; SRB

                                            ;;tmp_rain=(reflectivZ/aRRc)^(1/bRRc)                                                           ;; SRB
                                            ;;rain_R11Core[tmp_col,tmp_row,4]=rain_R11Core[tmp_col,tmp_row,4]+tmp_rain                      ;; SRB
                                            ;;nRai_R11Core[tmp_col,tmp_row,4]=nRai_R11Core[tmp_col,tmp_row,4]+1l                            ;; SRB
                                         endif else begin
                                            if raintypeFULL[col,fil] eq STRA then begin ;;stratiform rain
                                               RTo_stra_NSR_SH=RTo_stra_NSR_SH+nearSrfR_Org                  &  ctaStra_NSR_SH=ctaStra_NSR_SH+1l
                                               ;;RTo_stra_R11_SH=RTo_stra_R11_SH+(reflectivZ/aRRs)^(1/bRRs)    &  ctaStra_R11_SH=ctaStra_R11_SH+1l      ;; SRB
                                               
                                               ;;tmp_rain=(reflectivZ/aRRs)^(1/bRRs)                                                        ;; SRB
                                               ;;rain_R11Core[tmp_col,tmp_row,4]=rain_R11Core[tmp_col,tmp_row,4]+tmp_rain                   ;; SRB
                                               ;;nRai_R11Core[tmp_col,tmp_row,4]=nRai_R11Core[tmp_col,tmp_row,4]+1l                         ;; SRB
                                            endif else begin
                                               if raintypeFULL[col,fil] ge OTHER then begin ;;other rain
                                                  RTo_othe_NSR_SH=RTo_othe_NSR_SH+nearSrfR_Org               &  ctaOthe_NSR_SH=ctaOthe_NSR_SH+1l
                                                  ;;RTo_othe_R11_SH=RTo_othe_R11_SH+(reflectivZ/aRR)^(1/bRR)   &  ctaOthe_R11_SH=ctaOthe_R11_SH+1l      ;; SRB
                                        
                                                  ;;tmp_rain=(reflectivZ/aRR)^(1/bRR)                                                       ;; SRB
                                                  ;;rain_R11Core[tmp_col,tmp_row,4]=rain_R11Core[tmp_col,tmp_row,4]+tmp_rain                ;; SRB
                                                  ;;nRai_R11Core[tmp_col,tmp_row,4]=nRai_R11Core[tmp_col,tmp_row,4]+1l                      ;; SRB
                                               endif else begin
                                                  if raintype[col,fil] eq raintype_noRainValue then begin ;;No rain
                                                     RTo_noRa_NSR_SH=RTo_noRa_NSR_SH+0.   &   ctaNoRa_NSR_SH=ctaNoRa_NSR_SH+1l
                                                     ;;RTo_noRa_R11_SH=RTo_noRa_R11_SH+0.   &   ctaNoRa_R11_SH=ctaNoRa_R11_SH+1l                        ;; SRB
                                                     
                                                     ;;rain_R11Core[tmp_col,tmp_row,4]=rain_R11Core[tmp_col,tmp_row,4]+0                    ;; SRB
                                                     ;;nRai_R11Core[tmp_col,tmp_row,4]=nRai_R11Core[tmp_col,tmp_row,4]+1l                   ;; SRB
                                                  endif
                                               endelse
                                            endelse
                                         endelse
                                      endif    ;; end for missing values....
                                   endif       ;; if ctaC ne 0 and ctaR ne 0
                                ;;endif       ;; end for pixels with height > 2.5km                                                      ;; SRB
                             endif          ;;checking if there is data in the vertical
                          endif             ;;checking if there is data inside the boundaries of full storm  
                       endfor               ;;for going thry different pixelsum (indivdual pixels within storm)

                       ;;************************************************************************************************************    ;; SRB
                       ;;Romatschke and Houze 2011
                       ;;total_RainAll_R11=RTo_stra_R11_SH+RTo_conv_R11_SH+RTo_othe_R11_SH+RTo_noRa_R11_SH
                       ;;total_RainCSs_R11=RTo_stra_R11_SH+RTo_conv_R11_SH
                       ;;m_RainAll_R11=total_RainAll_R11/float(ctaStra_R11_SH+ctaConv_R11_SH+ctaOthe_R11_SH+ctaNoRa_R11_SH)  ;;mean rainfall for all pixels
                       ;;if ctaStra_R11_SH ne 0 then m_RainStrt_R11=RTo_stra_R11_SH/ctaStra_R11_SH else m_RainStrt_R11=0. ;;mean stratiform rain
                       ;;if ctaConv_R11_SH ne 0 then m_RainConv_R11=RTo_conv_R11_SH/ctaConv_R11_SH else m_RainConv_R11=0. ;;mean convective rain
                       ;;if ctaStra_R11_SH ne 0 or ctaConv_R11_SH ne 0 then $
                       ;;   m_RainCSs_R11=total_RainCSs_R11/float(ctaStra_R11_SH+ctaConv_R11_SH) else m_RainCSs_R11=-999.
                       ;;;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       ;;size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       ;;vol_Rain_All_R11=total_RainAll_R11*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Str_R11=RTo_stra_R11_SH*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;vol_Rain_Con_R11=RTo_conv_R11_SH*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       ;;tmp_statsRain_R11_SH=[m_RainAll_R11,m_RainStrt_R11,m_RainConv_R11,vol_Rain_All_R11,vol_Rain_Str_R11,vol_Rain_Con_R11]

                       ;;************************************************************************************************************
                       ;;Near Surface Rain
                       total_RainAll_NSR=RTo_stra_NSR_SH+RTo_conv_NSR_SH+RTo_othe_NSR_SH+RTo_noRa_NSR_SH
                       total_RainCSs_NSR=RTo_stra_NSR_SH+RTo_conv_NSR_SH
                       m_RainAll_NSR=total_RainAll_NSR/float(ctaStra_NSR_SH+ctaConv_NSR_SH+ctaOthe_NSR_SH+ctaNoRa_NSR_SH)  ;;mean rainfall for all pixels
                       if ctaStra_NSR_SH ne 0 then m_RainStrt_NSR=RTo_stra_NSR_SH/ctaStra_NSR_SH else m_RainStrt_NSR=0. ;;mean stratiform rain
                       if ctaConv_NSR_SH ne 0 then m_RainConv_NSR=RTo_conv_NSR_SH/ctaConv_NSR_SH else m_RainConv_NSR=0. ;;mean convective rain
                       if ctaStra_NSR_SH ne 0 or ctaConv_NSR_SH ne 0 then $
                          m_RainCSs_NSR=total_RainCSs_NSR/float(ctaStra_NSR_SH+ctaConv_NSR_SH) else m_RainCSs_NSR=-999.
                       ;;here calculate volumetric values of rain ;;volume in [1e6*kg/s]
                       size_pixels=deg2km(pixDeg,cen_lon,cen_lat)
                       vol_Rain_All_NSR=total_RainAll_NSR*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Str_NSR=RTo_stra_NSR_SH*(size_pixels[0]*size_pixels[1])/secsPerHr
                       vol_Rain_Con_NSR=RTo_conv_NSR_SH*(size_pixels[0]*size_pixels[1])/secsPerHr
                       ;;mean rain, mean strat, mean convec, vol all, vol strat, vol conv [1e6*kg/s]
                       tmp_statsRain_NSR_SH=[m_RainAll_NSR,m_RainStrt_NSR,m_RainConv_NSR,vol_Rain_All_NSR,vol_Rain_Str_NSR,vol_Rain_Con_NSR]

                       ;;*********************************************************************************************
                       ;;Here I calculate the CFAD count !!! for Full Storm!!! (only for one single Full storm)
                       if area ne areaIsFULL and cen_lon ne lonCIsFULL and cen_lat ne latCIsFull then begin
                          refl_SingleStorm=fltarr(nlonsFull,nlatsFull,nlevels)
                          refl_SingleStorm[*,*,*]=refl_3D_fillValue
                          refl_SingleStorm[w_idF]=refl_3D_FULL[w_idF]
                          if cta_Full ne npix_str[donde_shallow[ss]] then stop  ;; just to check! because this is in 3D!
                   
                          ;;here count reflectivity for each pixel that compose the storm into a matrix of CFAD
                          for i=0l,cta_Full-1l do begin  
                             col_R=where(refl_CFAD eq round(refl_3D_FULL[w_idF[i]]),ctaZ)   
                             row_H=where(alts_CFAD eq hgts_3D_FULL[w_idF[i]],ctaH)   ;;here locate the height of the pixel
                             if ctaH ne 0 and ctaZ ne 0 then CFAD_Full[col_R,row_H,4]=CFAD_Full[col_R,row_H,4]+1l   
                          endfor
                       endif

                       ;;*********************************************************************************************
                       ;;Here I calculate the CFAD count !!! for Convective component!!!
                       refl_SingleStorm=fltarr(nlonsFull,nlatsFull,nlevels)
                       refl_SingleStorm[*,*,*]=refl_3D_fillValue
                       refl_SingleStorm[w_id_SH]=refl_3D_FULL[w_id_SH]
                       if cta1_SH ne npix_SH[donde_Shallow2[ssSH]] then stop  ;; just to check! because this is in 3D!
                
                       ;;here count reflectivity for each pixel that compose the storm into a matrix of CFAD
                       for i=0l,cta1_SH-1 do begin  
                          col_R=where(refl_CFAD eq round(refl_3D_FULL[w_id_SH[i]]),ctaZ)   
                          row_H=where(alts_CFAD eq hgts_3D_FULL[w_id_SH[i]],ctaH)   ;;here locate the height of the pixel
                          if ctaH ne 0 and ctaZ ne 0 then CFAD_Core[col_R,row_H,4]=CFAD_Core[col_R,row_H,4]+1l
                       endfor

                       areaIsFULL=area & lonCIsFULL=cen_lon & latCIsFull=cen_lat   ;;return area counter to avoid double count of fullStorm

                       ;;info_SH=[info_SH,orbit+'.'+datetime+'.'+strtrim(string(s_id_SH),2)]  ;;MASK MOD
                       ;;info_SH=[info_SH,orbit+'.'+datetime+'.'+strtrim(string(num_shi),2)]  ;;MASK MOD
                       info_SH=[info_SH,orbit+'.'+datetime+'.0']  ;;MASK MOD
                       shape_Core_SH=[[shape_Core_SH],[tmp_shapeSH]]
                       shape_Full_SH=[[shape_Full_SH],[tmp_shape]]
                       rain_Core_SH=[[rain_Core_SH],[rain_momentSH]]
                       rain_Full_SH=[[rain_Full_SH],[rain_moment]]
                       rainTypeCore_SH=[[rainTypeCore_SH],[statsRain_SH]]
                       rainTypeFull_SH=[[rainTypeFull_SH],[statsRain]]
                
                       rainCore_SH_NSR=[[rainCore_SH_NSR],[tmp_statsRain_NSR_SH]]
                       rainFull_SH_NSR=[[rainFull_SH_NSR],[tmp_statsRain_NSR]]
                       ;;rainCore_SH_R11=[[rainCore_SH_R11],[tmp_statsRain_R11_SH]]                                                  ;; SRB
                       ;;rainFull_SH_R11=[[rainFull_SH_R11],[tmp_statsRain_R11]]                                                     ;; SRB

                       undefine,dim_lonSH
                       undefine,dim_latSH
                       undefine,hgt_sumSH
                       undefine,grid_sum
                       undefine,donde
                       undefine,pixelsum
                       undefine,lon_sum
                       undefine,lat_sum
                       undefine,hgt_sum
                       undefine,tmp_shapeSH
                       undefine,tmp_shape
                       undefine,SrfRain_FULL
                       undefine,raintypeFULL
                       undefine,refl_3D_FULL
                       undefine,hgts_3D_FULL
                       undefine,grid_storm_FULL
                       undefine,rainSH
                       undefine,rain
                       undefine,rain_nomiss
                       undefine,stratconv
                       undefine,strats
                       undefine,convec
                       undefine,others
                       undefine,noRain
                       undefine,missin
                       undefine,refl_SingleStorm
                       undefine,col_R
                       undefine,row_H
                    endif      ;;endif found a storm cluster with contiguous convective pixels within theresholds 
                    undefine,s_id_SH
                    undefine,w_id_SH
                    undefine,singlestormgrid_SH
                    undefine,grid_sum_SH_SH
                    undefine,dondeSH_SH
                    undefine,lonSH
                    undefine,latSH
                 endfor      ;;endfor loop through analyzed storms clusters that maybe are deep-wide convective
                 undefine,id_SH
                 undefine,npix_SH
                 undefine,grid_SH
                 undefine,searchNaN_SH
                 undefine,donde_Shallow2
              endif      ;;endif for convective areas that could be matching the theresholds
              undefine,lonSH
              undefine,latSH
              undefine,dondeSH
              undefine,grid_sumSH
              undefine,singlestormgridShallow
              undefine,lonsFull_sub
              undefine,latsFull_sub
              undefine,singlestormgrid_Full
              undefine,d_latsFull
              undefine,d_lonsFull
              undefine,total_latFull
              undefine,total_lonFull
              undefine,s_idF
              undefine,w_idF
           endfor      ;;endfor loop thru full storm volumens greater than 2l
        endif       ;;end if the orbit does have shallow isolated pixels

        undefine,grid_storm
        undefine,npix_str
        undefine,id_storm
        undefine,searchNaN
     endif    ;;end if of pixels with reflectivity greater than 0 (its raining) ;;***********
     undefine,refl_Rain
     undefine,hgts_3D
     undefine,rain_type3D
     undefine,d_strt
     undefine,d_conv
     undefine,d_Shal
     undefine,d_ShIs
     undefine,d_ShNi
     undefine,d_othe
     undefine,d_noRa
     undefine,lons
     undefine,lats
     undefine,hgts
     undefine,refl_3D
     undefine,SrfRain
     undefine,raintype
     undefine,lons3D
     undefine,lats3D

     ;; Write mask info to input file and deallocate mask arrays
     if makeCoreMasks then begin
        bsr_id = ncdf_varid(ncid,'bsr_mask_str')
        if bsr_id ne -1 then ncdf_varput,ncid,bsr_id,bsr_mask
        dcc_id = ncdf_varid(ncid,'dcc_mask_str')
        if dcc_id ne -1 then ncdf_varput,ncid,dcc_id,dcc_mask
        dwc_id = ncdf_varid(ncid,'dwc_mask_str')
        if dwc_id ne -1 then ncdf_varput,ncid,dwc_id,dwc_mask
        wcc_id = ncdf_varid(ncid,'wcc_mask_str')
        if wcc_id ne -1 then ncdf_varput,ncid,wcc_id,wcc_mask
        ;;shi_id = ncdf_varid(ncid,'shi_mask_str')
        ;;if shi_id ne -1 then ncdf_varput,ncid,shi_id,shi_mask
        storm_id = ncdf_varid(ncid,'storm_mask_str')
        if storm_id ne -1 then ncdf_varput,ncid,storm_id,storm_mask

        undefine,bsr_mask
        undefine,dcc_mask
        undefine,dwc_mask
        undefine,wcc_mask
        ;;undefine,shi_mask
        undefine,storm_mask
     endif 

     ;; close input file
     ncdf_close,ncid
     
  endfor     ;;end for multiples file orbits within the directory (year-month)

  ;;*********************************************************************************************************************
  ;;Here it will save the individual monthly files with detailed info of each core and full storm - check for directories
  if file_test(path_out+'monthly_class_v11s',/directory) eq 0 then file_mkdir,path_out+'monthly_class_v11s'
  if file_test(path_out+'monthly_class_v11s/'+meses,/directory) eq 0 then file_mkdir,path_out+'monthly_class_v11s/'+meses
  
  ;;Save the Deep Convective core data  ;;****************************************************************************
  openw,1,path_out+'monthly_class_v11s/'+meses+'/Deep_Convective_'+type+'_v11s.dat'
  for i=1l,n_elements(info_DC)-1 do $
     printf,1,format='(9f10.2,7f10.2,6f11.4,9f10.2,7f10.2,6f11.4)',$
            shape_Core_DC[*,i],rain_Core_DC[*,i],rainTypeCore_DC[*,i],shape_Full_DC[*,i],rain_Full_DC[*,i],rainTypeFull_DC[*,i]
  close,1
  openw,1,path_out+'monthly_class_v11s/'+meses+'/Deep_Convective_'+type+'_v11s.info'
  ;; MASK MOD - change format statment to max len of a record in info_DC
  ;;for i=1l,n_elements(info_DC)-1 do printf,1,format='(a25)',info_DC[i]
  for i=1l,n_elements(info_DC)-1 do begin
     alen = '(a'+strtrim(string(strlen(info_DC[i])),2)+')'
     printf,1,format=alen,info_DC[i]
  endfor 
  close,1

  ;;Save the Wide Convective core data  ;;****************************************************************************
  openw,1,path_out+'monthly_class_v11s/'+meses+'/Wide_Convective_'+type+'_v11s.dat'
  for i=1l,n_elements(info_WC)-1 do $
     printf,1,format='(9f10.2,7f10.2,6f11.4,9f10.2,7f10.2,6f11.4)',$
            shape_Core_WC[*,i],rain_Core_WC[*,i],rainTypeCore_WC[*,i],shape_Full_WC[*,i],rain_Full_WC[*,i],rainTypeFull_WC[*,i]
  close,1
  openw,1,path_out+'monthly_class_v11s/'+meses+'/Wide_Convective_'+type+'_v11s.info'
  ;; MASK MOD - change format statment to max len of a record in info_WC
  ;;for i=1l,n_elements(info_WC)-1 do printf,1,format='(a25)',info_WC[i]
  for i=1l,n_elements(info_WC)-1 do begin
     alen = '(a'+strtrim(string(strlen(info_WC[i])),2)+')'
     printf,1,format=alen,info_WC[i]
  endfor 
  close,1

  ;;Save the Deep and Wide Convective core data  ;;************************************************************************
  openw,1,path_out+'monthly_class_v11s/'+meses+'/DeepWide_Convective_'+type+'_v11s.dat'
  for i=1l,n_elements(info_DW)-1 do $
     printf,1,format='(9f10.2,7f10.2,6f11.4,9f10.2,7f10.2,6f11.4)',$
            shape_Core_DW[*,i],rain_Core_DW[*,i],rainTypeCore_DW[*,i],shape_Full_DW[*,i],rain_Full_DW[*,i],rainTypeFull_DW[*,i]
  close,1
  openw,1,path_out+'monthly_class_v11s/'+meses+'/DeepWide_Convective_'+type+'_v11s.info'
  ;; MASK MOD - change format statment to max len of a record in info_DW
  ;;for i=1l,n_elements(info_DW)-1 do printf,1,format='(a25)',info_DW[i]
  for i=1l,n_elements(info_DW)-1 do begin
     alen = '(a'+strtrim(string(strlen(info_DW[i])),2)+')'
     printf,1,format=alen,info_DW[i]
  endfor 
  close,1

  ;;Save the Broad Stratiform Regions data  ;;****************************************************************************
  openw,1,path_out+'monthly_class_v11s/'+meses+'/BroadStratiform_'+type+'_v11s.dat'
  for i=1l,n_elements(info_BS)-1 do $
     printf,1,format='(9f10.2,7f10.2,6f11.4,9f10.2,7f10.2,6f11.4)',$
            shape_Core_BS[*,i],rain_Core_BS[*,i],rainTypeCore_BS[*,i],shape_Full_BS[*,i],rain_Full_BS[*,i],rainTypeFull_BS[*,i]
  close,1
  openw,1,path_out+'monthly_class_v11s/'+meses+'/BroadStratiform_'+type+'_v11s.info'
  ;; MASK MOD - change format statment to max len of a record in info_BS
  ;;for i=1l,n_elements(info_BS)-1 do printf,1,format='(a25)',info_BS[i]
  for i=1l,n_elements(info_BS)-1 do begin
     alen = '(a'+strtrim(string(strlen(info_BS[i])),2)+')'
     printf,1,format=alen,info_BS[i]
  endfor 
  close,1

  ;;Save the Shallow Isolated events  ;;****************************************************************************
  openw,1,path_out+'monthly_class_v11s/'+meses+'/ShallowIsolated_'+type+'_v11s.dat'
  for i=1l,n_elements(info_SH)-1 do $
     printf,1,format='(9f10.2,7f10.2,6f11.4,9f10.2,7f10.2,6f11.4)',$
            shape_Core_SH[*,i],rain_Core_SH[*,i],rainTypeCore_SH[*,i],shape_Full_SH[*,i],rain_Full_SH[*,i],rainTypeFull_SH[*,i]
  close,1
  openw,1,path_out+'monthly_class_v11s/'+meses+'/ShallowIsolated_'+type+'_v11s.info'
  ;; MASK MOD - change format statment to max len of a record in info_SH
  for i=1l,n_elements(info_SH)-1 do printf,1,format='(a24)',info_SH[i]
  ;;for i=1l,n_elements(info_SH)-1 do begin
     ;;alen = '(a'+strtrim(string(strlen(info_SH[i])),2)+')'
     ;;printf,1,format=alen,info_SH[i]
  ;;endfor 
  close,1
  ;;*********************************************************************************************************************
     
  ;;*********************************************************************************************************************
  ;;Near Surface Rain  ;;here it saves: cen_lon,cen_lat,area_storm,max_height - check for directories
  if file_test(path_out+'stats_class_v11s',/directory) eq 0 then file_mkdir,path_out+'stats_class_v11s'
  if file_test(path_out+'stats_class_v11s/NearSurfRain',/directory) eq 0 then file_mkdir,path_out+'stats_class_v11s/NearSurfRain'
  if file_test(path_out+'stats_class_v11s/NearSurfRain/'+meses,/directory) eq 0 then  $
     file_mkdir,path_out+'stats_class_v11s/NearSurfRain/'+meses
  
  ;;Save the Deep Convective core data  ;;************************************************************
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Deep_Convec_'+type+'_v11s.stats'
  printf,34,n_elements(info_DC)-1
  for i=1l,n_elements(info_DC)-1 do $
     printf,34,format='(4f12.3,4f12.3,a27)',$
            shape_Core_DC[0:3,i],shape_Full_DC[0:3,i],info_DC[i]
  close,34

  ;;here it saves: mean_rain,mean_strat,mean_convec,vol_all,vol_strat,vol_conv  [1e6*kg/s] 
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Deep_Convec_'+type+'_v11s.rain'
  printf,34,n_elements(info_DC)-1
  for i=1,n_elements(info_DC)-1 do printf,34,format='(24f12.2,a27)',$
                                          rainCore_DC_NSR[*,i],rainFull_DC_NSR[*,i],rainTypeCore_DC[*,i],rainTypeFull_DC[*,i],info_DC[i]
  close,34

  ;;Save the Wide Convective core data   ;;****************************************************************************
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Wide_Convec_'+type+'_v11s.stats'
  printf,34,n_elements(info_WC)-1
  for i=1l,n_elements(info_WC)-1 do $
     printf,34,format='(4f12.3,4f12.3,a27)',$
            shape_Core_WC[0:3,i],shape_Full_WC[0:3,i],info_WC[i]
  close,34

  ;;here it saves: mean_rain,mean_strat,mean_convec,vol_all,vol_strat,vol_conv  [1e6*kg/s] 
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Wide_Convec_'+type+'_v11s.rain'
  printf,34,n_elements(info_WC)-1
  for i=1,n_elements(info_WC)-1 do printf,34,format='(24f12.2,a27)',$
                                          rainCore_WC_NSR[*,i],rainFull_WC_NSR[*,i],rainTypeCore_WC[*,i],rainTypeFull_WC[*,i],info_WC[i]
  close,34

  ;;Save the Deep and Wide Convective core data  ;;***********************************************************************
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/DeepWide_Convec_'+type+'_v11s.stats'
  printf,34,n_elements(info_DW)-1
  for i=1l,n_elements(info_DW)-1 do $
     printf,34,format='(4f12.3,4f12.3,a27)',$
            shape_Core_DW[0:3,i],shape_Full_DW[0:3,i],info_DW[i]
  close,34

  ;;here it saves: mean_rain,mean_strat,mean_convec,vol_all,vol_strat,vol_conv  [1e6*kg/s] 
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/DeepWide_Convec_'+type+'_v11s.rain'
  printf,34,n_elements(info_DW)-1
  for i=1,n_elements(info_DW)-1 do printf,34,format='(24f12.2,a27)',$
                                          rainCore_DW_NSR[*,i],rainFull_DW_NSR[*,i],rainTypeCore_DW[*,i],rainTypeFull_DW[*,i],info_DW[i]
  close,34

  ;;Save the Broad Stratiform Regions data  ;;****************************************************************************
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Broad_Strat_'+type+'_v11s.stats'
  printf,34,n_elements(info_BS)-1
  for i=1l,n_elements(info_BS)-1 do $
     printf,34,format='(4f12.3,4f12.3,a27)',$
            shape_Core_BS[0:3,i],shape_Full_BS[0:3,i],info_BS[i]
  close,34

  ;;here it saves: mean_rain,mean_strat,mean_convec,vol_all,vol_strat,vol_conv  [1e6*kg/s] 
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Broad_Strat_'+type+'_v11s.rain'
  printf,34,n_elements(info_BS)-1
  for i=1,n_elements(info_BS)-1 do printf,34,format='(24f12.2,a27)',$
                                          rainCore_BS_NSR[*,i],rainFull_BS_NSR[*,i],rainTypeCore_BS[*,i],rainTypeFull_BS[*,i],info_BS[i]
  close,34

  ;;Save the Shallow Isolated data  ;;****************************************************************************
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Shallow_Isol_'+type+'_v11s.stats'
  printf,34,n_elements(info_SH)-1
  for i=1l,n_elements(info_SH)-1 do $
     printf,34,format='(4f12.3,4f12.3,a27)',$
            shape_Core_SH[0:3,i],shape_Full_SH[0:3,i],info_SH[i]
  close,34

  ;;here it saves: mean_rain,mean_strat,mean_convec,vol_all,vol_strat,vol_conv  [1e6*kg/s] 
  openw,34,path_out+'stats_class_v11s/NearSurfRain/'+meses+'/Shallow_Isol_'+type+'_v11s.rain'
  printf,34,n_elements(info_SH)-1
  for i=1,n_elements(info_SH)-1 do printf,34,format='(24f12.2,a27)',$
                                          rainCore_SH_NSR[*,i],rainFull_SH_NSR[*,i],rainTypeCore_SH[*,i],rainTypeFull_SH[*,i],info_SH[i]
  close,34
  ;;*********************************************************************************************************************

  ;;*********************************************************************************************************************
  ;;save the Rain accumulation Matrices for Near Surface Rain!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name=path_out+'stats_class_v11s/NearSurfRain/'+meses+'/infoRainNSR_Stra_EchoCores_'+type+'_v11s'
  id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ncdf_control,id,/fill                                  ;; Fill the file with default values
  lnC_id = ncdf_dimdef(id,'lons_c',nlonsC)
  ltC_id = ncdf_dimdef(id,'lats_c',nlatsC)
  lnF_id = ncdf_dimdef(id,'lons_f',nlonsF)
  ltF_id = ncdf_dimdef(id,'lats_f',nlatsF)

  freqF_id=ncdf_vardef(id,'freq_Full',[lnC_id,ltC_id],/long) ;; Define variables:
  ncdf_attput,id,freqF_id,'long_name','pixel count for full storm',/char
  ncdf_attput,id,freqF_id,'missing_value','-9999.',/char

  freqC_id=ncdf_vardef(id,'freq_core',[lnC_id,ltC_id],/long) ;; Define variables:
  ncdf_attput,id,freqC_id,'long_name','pixel count for core storm',/char
  ncdf_attput,id,freqC_id,'missing_value','-9999.',/char
  ;;**    
  rainF_id_BS=ncdf_vardef(id,'rain_Full_BS',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rainF_id_BS,'long_name','Accumulated rain for full storm in BS',/char
  ncdf_attput,id,rainF_id_BS,'missing_value','-9999.',/char

  nR_F_id_BS=ncdf_vardef(id,'nRain_Full_BS',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_F_id_BS,'long_name','number of elements of rain for full storm in BS',/char
  ncdf_attput,id,nR_F_id_BS,'missing_value','-9999.',/char

  rain_id_BS=ncdf_vardef(id,'rain_Core_BS',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rain_id_BS,'long_name','Accumulated rain for BS cores',/char
  ncdf_attput,id,rain_id_BS,'missing_value','-9999.',/char

  nR_id_BS=ncdf_vardef(id,'nRain_Core_BS',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_id_BS,'long_name','number of elements of rain for BS cores',/char
  ncdf_attput,id,nR_id_BS,'missing_value','-9999.',/char

  lonID_C = ncdf_vardef(id,'lonsC',[lnC_id],/float)
  ncdf_attput,id,lonID_C,'units','degrees',/char
  ncdf_attput,id,lonID_C,'long_name','longitudes in coarse grid',/char
  latID_C = ncdf_vardef(id,'latsC',[ltC_id],/float)
  ncdf_attput,id,latID_C,'units','degrees',/char
  ncdf_attput,id,latID_C,'long_name','latitudes in coarse grid',/char

  lonID_F = ncdf_vardef(id,'lonsF',[lnF_id],/double)
  ncdf_attput,id,lonID_F,'units','degrees',/char
  ncdf_attput,id,lonID_F,'long_name','longitudes in fine grid',/char
  latID_F = ncdf_vardef(id,'latsF',[ltF_id],/double)
  ncdf_attput,id,latID_F,'units','degrees',/char
  ncdf_attput,id,latID_F,'long_name','latitudes in fine grid',/char

  ncdf_attput,id,/global,'Title','Information of Broad Stratiform classified storms'
  ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ncdf_attput,id,/global,'rain_text1','rain calculated using Near Surface Rain'
  ncdf_control,id,/endef                    ;; Put file in data mode:

  ncdf_varput,id,freqF_id,freq_Full[*,*,3]                         ;;Write the data into file
  ncdf_varput,id,freqC_id,freq_Core[*,*,3]

  ncdf_varput,id,rainF_id_BS,rain_NSRFull[*,*,3]
  ncdf_varput,id,nR_F_id_BS,nRai_NSRFull[*,*,3]
  ncdf_varput,id,rain_id_BS,rain_NSRCore[*,*,3]
  ncdf_varput,id,nR_id_BS,nRai_NSRCore[*,*,3]
  
  ncdf_varput,id,lonID_C,lonsC
  ncdf_varput,id,latID_C,latsC
  ncdf_varput,id,lonID_F,lonsF
  ncdf_varput,id,latID_F,latsF

  ncdf_close,id
  spawn,'zip -mjq '+name+'.zip '+name+'.nc'

  ;;save the accumulated rain matrices in Convective Cores computed using NEar surface rain
  name=path_out+'stats_class_v11s/NearSurfRain/'+meses+'/infoRainNSR_Conv_EchoCores_'+type+'_v11s'
  id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ncdf_control,id,/fill                                  ;; Fill the file with default values
  lnC_id = ncdf_dimdef(id,'lons_c',nlonsC)
  ltC_id = ncdf_dimdef(id,'lats_c',nlatsC)
  lnF_id = ncdf_dimdef(id,'lons_f',nlonsF)
  ltF_id = ncdf_dimdef(id,'lats_f',nlatsF)
  var_id = ncdf_dimdef(id,'echo_type',3)

  freqF_id=ncdf_vardef(id,'freq_Full',[lnC_id,ltC_id,var_id],/long) ;; Define variables:
  ncdf_attput,id,freqF_id,'long_name','Event pixel count for Entire storm',/char
  ncdf_attput,id,freqF_id,'missing_value','-9999.',/char

  freqC_id=ncdf_vardef(id,'freq_core',[lnC_id,ltC_id,var_id],/long) ;; Define variables:
  ncdf_attput,id,freqC_id,'long_name','Event pixel count for cores',/char
  ncdf_attput,id,freqC_id,'missing_value','-9999.',/char
  ;;** 
  rainF_id_DC=ncdf_vardef(id,'rain_Full_DC',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rainF_id_DC,'long_name','Accumulated rain for Entire storm in DC',/char
  ncdf_attput,id,rainF_id_DC,'missing_value','-9999.',/char

  nR_F_id_DC=ncdf_vardef(id,'nRain_Full_DC',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_F_id_DC,'long_name','number of elements of rain for Entire storm in DC',/char
  ncdf_attput,id,nR_F_id_DC,'missing_value','-9999.',/char

  rain_id_DC=ncdf_vardef(id,'rain_Core_DC',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rain_id_DC,'long_name','Accumulated rain for DC cores',/char
  ncdf_attput,id,rain_id_DC,'missing_value','-9999.',/char

  nR_id_DC=ncdf_vardef(id,'nRain_Core_DC',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_id_DC,'long_name','number of elements of rain for DC cores',/char
  ncdf_attput,id,nR_id_DC,'missing_value','-9999.',/char
  ;;** 
  rainF_id_WC=ncdf_vardef(id,'rain_Full_WC',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rainF_id_WC,'long_name','Accumulated rain for Entire storm in WC',/char
  ncdf_attput,id,rainF_id_WC,'missing_value','-9999.',/char

  nR_F_id_WC=ncdf_vardef(id,'nRain_Full_WC',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_F_id_WC,'long_name','number of elements of rain for Entire storm in WC',/char
  ncdf_attput,id,nR_F_id_WC,'missing_value','-9999.',/char

  rain_id_WC=ncdf_vardef(id,'rain_Core_WC',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rain_id_WC,'long_name','Accumulated rain for WC cores',/char
  ncdf_attput,id,rain_id_WC,'missing_value','-9999.',/char

  nR_id_WC=ncdf_vardef(id,'nRain_Core_WC',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_id_WC,'long_name','number of elements of rain for WC cores',/char
  ncdf_attput,id,nR_id_WC,'missing_value','-9999.',/char
  ;;** 
  rainF_id_DW=ncdf_vardef(id,'rain_Full_DW',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rainF_id_DW,'long_name','Accumulated rain for Entire storm in DW',/char
  ncdf_attput,id,rainF_id_DW,'missing_value','-9999.',/char

  nR_F_id_DW=ncdf_vardef(id,'nRain_Full_DW',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_F_id_DW,'long_name','number of elements of rain for Entire storm in DW',/char
  ncdf_attput,id,nR_F_id_DW,'missing_value','-9999.',/char

  rain_id_DW=ncdf_vardef(id,'rain_Core_DW',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rain_id_DW,'long_name','Accumulated rain for DW cores',/char
  ncdf_attput,id,rain_id_DW,'missing_value','-9999.',/char

  nR_id_DW=ncdf_vardef(id,'nRain_Core_DW',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_id_DW,'long_name','number of elements of rain for DW cores',/char
  ncdf_attput,id,nR_id_DW,'missing_value','-9999.',/char
  ;;**
  lonID_C = ncdf_vardef(id,'lonsC',[lnC_id],/float)
  ncdf_attput,id,lonID_C,'units','degrees',/char
  ncdf_attput,id,lonID_C,'long_name','longitudes in coarse grid',/char
  latID_C = ncdf_vardef(id,'latsC',[ltC_id],/float)
  ncdf_attput,id,latID_C,'units','degrees',/char
  ncdf_attput,id,latID_C,'long_name','latitudes in coarse grid',/char

  lonID_F = ncdf_vardef(id,'lonsF',[lnF_id],/double)
  ncdf_attput,id,lonID_F,'units','degrees',/char
  ncdf_attput,id,lonID_F,'long_name','longitudes in fine grid',/char
  latID_F = ncdf_vardef(id,'latsF',[ltF_id],/double)
  ncdf_attput,id,latID_F,'units','degrees',/char
  ncdf_attput,id,latID_F,'long_name','latitudes in fine grid',/char

  ncdf_attput,id,/global,'Title','Information of Convective (DCC-WCC-DWCC) classified storms'
  ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ncdf_attput,id,/global,'rain_text1','rain calculated using Near Surface Rain'
  ncdf_control,id,/endef                    ;; Put file in data mode:

  ncdf_varput,id,freqF_id,freq_Full[*,*,[0,1,2]]                         ;;Write the data into file
  ncdf_varput,id,freqC_id,freq_Core[*,*,[0,1,2]]

  ncdf_varput,id,rainF_id_DC,rain_NSRFull[*,*,0]
  ncdf_varput,id,nR_F_id_DC,nRai_NSRFull[*,*,0]
  ncdf_varput,id,rain_id_DC,rain_NSRCore[*,*,0]
  ncdf_varput,id,nR_id_DC,nRai_NSRCore[*,*,0]

  ncdf_varput,id,rainF_id_WC,rain_NSRFull[*,*,1]
  ncdf_varput,id,nR_F_id_WC,nRai_NSRFull[*,*,1]
  ncdf_varput,id,rain_id_WC,rain_NSRCore[*,*,1]
  ncdf_varput,id,nR_id_WC,nRai_NSRCore[*,*,1]
  
  ncdf_varput,id,rainF_id_DW,rain_NSRFull[*,*,2]
  ncdf_varput,id,nR_F_id_DW,nRai_NSRFull[*,*,2]
  ncdf_varput,id,rain_id_DW,rain_NSRCore[*,*,2]
  ncdf_varput,id,nR_id_DW,nRai_NSRCore[*,*,2]

  ncdf_varput,id,lonID_C,lonsC
  ncdf_varput,id,latID_C,latsC
  ncdf_varput,id,lonID_F,lonsF
  ncdf_varput,id,latID_F,latsF
  
  ncdf_close,id
  spawn,'zip -mjq '+name+'.zip '+name+'.nc'
  ;;*********************************************************************************************************************

  ;;*********************************************************************************************************************
  ;;save the Rain accumulation Matrices for Shallow Isolated with Near Surface Rain!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name=path_out+'stats_class_v11s/NearSurfRain/'+meses+'/infoRainNSR_ShallowIsol_EchoCores_'+type+'_v11s'
  id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ncdf_control,id,/fill                                  ;; Fill the file with default values
  lnC_id = ncdf_dimdef(id,'lons_c',nlonsC)
  ltC_id = ncdf_dimdef(id,'lats_c',nlatsC)
  lnF_id = ncdf_dimdef(id,'lons_f',nlonsF)
  ltF_id = ncdf_dimdef(id,'lats_f',nlatsF)

  freqF_id=ncdf_vardef(id,'freq_Full',[lnC_id,ltC_id],/long) ;; Define variables:
  ncdf_attput,id,freqF_id,'long_name','pixel count for full storm',/char
  ncdf_attput,id,freqF_id,'missing_value','-9999.',/char

  freqC_id=ncdf_vardef(id,'freq_core',[lnC_id,ltC_id],/long) ;; Define variables:
  ncdf_attput,id,freqC_id,'long_name','pixel count for core storm',/char
  ncdf_attput,id,freqC_id,'missing_value','-9999.',/char
  ;;**    
  rainF_id_SH=ncdf_vardef(id,'rain_Full_SH',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rainF_id_SH,'long_name','Accumulated rain for full storm in SH',/char
  ncdf_attput,id,rainF_id_SH,'missing_value','-9999.',/char

  nR_F_id_SH=ncdf_vardef(id,'nRain_Full_SH',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_F_id_SH,'long_name','number of elements of rain for full storm in SH',/char
  ncdf_attput,id,nR_F_id_SH,'missing_value','-9999.',/char

  rain_id_SH=ncdf_vardef(id,'rain_Core_SH',[lnF_id,ltF_id],/float) ;; Define variables:
  ncdf_attput,id,rain_id_SH,'long_name','Accumulated rain for SH cores',/char
  ncdf_attput,id,rain_id_SH,'missing_value','-9999.',/char

  nR_id_SH=ncdf_vardef(id,'nRain_Core_SH',[lnF_id,ltF_id],/long) ;; Define variables:
  ncdf_attput,id,nR_id_SH,'long_name','number of elements of rain for SH cores',/char
  ncdf_attput,id,nR_id_SH,'missing_value','-9999.',/char

  lonID_C = ncdf_vardef(id,'lonsC',[lnC_id],/float)
  ncdf_attput,id,lonID_C,'units','degrees',/char
  ncdf_attput,id,lonID_C,'long_name','longitudes in coarse grid',/char
  latID_C = ncdf_vardef(id,'latsC',[ltC_id],/float)
  ncdf_attput,id,latID_C,'units','degrees',/char
  ncdf_attput,id,latID_C,'long_name','latitudes in coarse grid',/char

  lonID_F = ncdf_vardef(id,'lonsF',[lnF_id],/double)
  ncdf_attput,id,lonID_F,'units','degrees',/char
  ncdf_attput,id,lonID_F,'long_name','longitudes in fine grid',/char
  latID_F = ncdf_vardef(id,'latsF',[ltF_id],/double)
  ncdf_attput,id,latID_F,'units','degrees',/char
  ncdf_attput,id,latID_F,'long_name','latitudes in fine grid',/char

  ncdf_attput,id,/global,'Title','Information of Shallow Isolated cores'
  ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ncdf_attput,id,/global,'rain_text1','rain calculated using Near Surface Rain'
  ncdf_control,id,/endef                    ;; Put file in data mode:

  ncdf_varput,id,freqF_id,freq_Full[*,*,4]                         ;;Write the data into file
  ncdf_varput,id,freqC_id,freq_Core[*,*,4]

  ncdf_varput,id,rainF_id_SH,rain_NSRFull[*,*,4]
  ncdf_varput,id,nR_F_id_SH,nRai_NSRFull[*,*,4]
  ncdf_varput,id,rain_id_SH,rain_NSRCore[*,*,4]
  ncdf_varput,id,nR_id_SH,nRai_NSRCore[*,*,4]

  ncdf_varput,id,lonID_C,lonsC
  ncdf_varput,id,latID_C,latsC
  ncdf_varput,id,lonID_F,lonsF
  ncdf_varput,id,latID_F,latsF

  ncdf_close,id
  spawn,'zip -mjq '+name+'.zip '+name+'.nc'
  ;;*********************************************************************************************************************
  
  ;;*********************************************************************************************************************   ;; SRB
  ;;;;here it saves: cen_lon,cen_lat,area_storm,max_height
  ;;;;Romatschke and Houze 2011 ZR
  ;;if file_test(path_out+'stats_class_v11s/RomHouze_2011',/directory) eq 0 then file_mkdir,path_out+'stats_class_v11s/RomHouze_2011'
  ;;if file_test(path_out+'stats_class_v11s/RomHouze_2011/'+meses,/directory) eq 0 then  $
  ;;   file_mkdir,path_out+'stats_class_v11s/RomHouze_2011/'+meses
  
  ;;;;Save the Deep Convective core data  ;;*******************************************
  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Deep_Convec_'+type+'_v11s.stats'
  ;;printf,34,n_elements(info_DC)-1
  ;;for i=1l,n_elements(info_DC)-1 do $
  ;;   printf,34,format='(4f12.3,4f12.3,a27)',$
  ;;          shape_Core_DC[0:3,i],shape_Full_DC[0:3,i],info_DC[i]
  ;;close,34

  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Deep_Convec_'+type+'_v11s.rain'
  ;;printf,34,n_elements(info_DC)-1
  ;;for i=1,n_elements(info_DC)-1 do printf,34,format='(24f12.2,a27)',$
  ;;                                        rainCore_DC_R11[*,i],rainFull_DC_R11[*,i],rainTypeCore_DC[*,i],rainTypeFull_DC[*,i],info_DC[i]
  ;;close,34

  ;;;;Save the Wide Convective core data  ;;*******************************************
  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Wide_Convec_'+type+'_v11s.stats'
  ;;printf,34,n_elements(info_WC)-1
  ;;for i=1l,n_elements(info_WC)-1 do $
  ;;   printf,34,format='(4f12.3,4f12.3,a27)',$
  ;;          shape_Core_WC[0:3,i],shape_Full_WC[0:3,i],info_WC[i]
  ;;close,34

  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Wide_Convec_'+type+'_v11s.rain'
  ;;printf,34,n_elements(info_WC)-1
  ;;for i=1,n_elements(info_WC)-1 do printf,34,format='(24f12.2,a27)',$
  ;;                                        rainCore_WC_R11[*,i],rainFull_WC_R11[*,i],rainTypeCore_WC[*,i],rainTypeFull_WC[*,i],info_WC[i]
  ;;close,34

  ;;;;Save the Deep and Wide Convective core data  ;;*******************************************
  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/DeepWide_Convec_'+type+'_v11s.stats'
  ;;printf,34,n_elements(info_DW)-1
  ;;for i=1l,n_elements(info_DW)-1 do $
  ;;   printf,34,format='(4f12.3,4f12.3,a27)',$
  ;;          shape_Core_DW[0:3,i],shape_Full_DW[0:3,i],info_DW[i]
  ;;close,34

  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/DeepWide_Convec_'+type+'_v11s.rain'
  ;;printf,34,n_elements(info_DW)-1
  ;;for i=1,n_elements(info_DW)-1 do printf,34,format='(24f12.2,a27)',$
  ;;                                        rainCore_DW_R11[*,i],rainFull_DW_R11[*,i],rainTypeCore_DW[*,i],rainTypeFull_DW[*,i],info_DW[i]
  ;;close,34

  ;;;;Save the Broad Stratiform Regions data  ;;*******************************************
  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Broad_Strat_'+type+'_v11s.stats'
  ;;printf,34,n_elements(info_BS)-1
  ;;for i=1l,n_elements(info_BS)-1 do $
  ;;   printf,34,format='(4f12.3,4f12.3,a27)',$
  ;;          shape_Core_BS[0:3,i],shape_Full_BS[0:3,i],info_BS[i]
  ;;close,34

  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Broad_Strat_'+type+'_v11s.rain'
  ;;printf,34,n_elements(info_BS)-1
  ;;for i=1,n_elements(info_BS)-1 do printf,34,format='(24f12.2,a27)',$
  ;;                                        rainCore_BS_R11[*,i],rainFull_BS_R11[*,i],rainTypeCore_BS[*,i],rainTypeFull_BS[*,i],info_BS[i]
  ;;close,34

  ;;;;Save the Shallow Isolated data  ;;*******************************************
  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Shallow_Isol_'+type+'_v11s.stats'
  ;;printf,34,n_elements(info_SH)-1
  ;;for i=1l,n_elements(info_SH)-1 do $
  ;;   printf,34,format='(4f12.3,4f12.3,a27)',$
  ;;          shape_Core_SH[0:3,i],shape_Full_SH[0:3,i],info_SH[i]
  ;;close,34

  ;;openw,34,path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/Shallow_Isol_'+type+'_v11s.rain'
  ;;printf,34,n_elements(info_SH)-1
  ;;for i=1,n_elements(info_SH)-1 do printf,34,format='(24f12.2,a27)',$
  ;;                                        rainCore_SH_R11[*,i],rainFull_SH_R11[*,i],rainTypeCore_SH[*,i],rainTypeFull_SH[*,i],info_SH[i]
  ;;close,34
  ;;;;*********************************************************************************************************************

  ;;*********************************************************************************************************************            ;; SRB
  ;;;;save the Rain accumulation Matrices using Romatschke and Houze 2011!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;;name=path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/infoRainR11_Stra_EchoCores_'+type+'_v11s'
  ;;id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ;;ncdf_control,id,/fill                                  ;; Fill the file with default values
  ;;lnC_id = ncdf_dimdef(id,'lons_c',nlonsC)
  ;;ltC_id = ncdf_dimdef(id,'lats_c',nlatsC)
  ;;lnF_id = ncdf_dimdef(id,'lons_f',nlonsF)
  ;;ltF_id = ncdf_dimdef(id,'lats_f',nlatsF)

  ;;freqF_id=ncdf_vardef(id,'freq_Full',[lnC_id,ltC_id],/long) ;; Define variables:
  ;;ncdf_attput,id,freqF_id,'long_name','pixel count for full storm',/char
  ;;ncdf_attput,id,freqF_id,'missing_value','-9999.',/char

  ;;freqC_id=ncdf_vardef(id,'freq_core',[lnC_id,ltC_id],/long) ;; Define variables:
  ;;ncdf_attput,id,freqC_id,'long_name','pixel count for core storm',/char
  ;;ncdf_attput,id,freqC_id,'missing_value','-9999.',/char
  ;;;;**    
  ;;rainF_id_BS=ncdf_vardef(id,'rain_Full_BS',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rainF_id_BS,'long_name','Accumulated rain for full storm in BS',/char
  ;;ncdf_attput,id,rainF_id_BS,'missing_value','-9999.',/char

  ;;nR_F_id_BS=ncdf_vardef(id,'nRain_Full_BS',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_F_id_BS,'long_name','number of elements of rain for full storm in BS',/char
  ;;ncdf_attput,id,nR_F_id_BS,'missing_value','-9999.',/char

  ;;rain_id_BS=ncdf_vardef(id,'rain_Core_BS',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rain_id_BS,'long_name','Accumulated rain for BS cores',/char
  ;;ncdf_attput,id,rain_id_BS,'missing_value','-9999.',/char

  ;;nR_id_BS=ncdf_vardef(id,'nRain_Core_BS',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_id_BS,'long_name','number of elements of rain for BS cores',/char
  ;;ncdf_attput,id,nR_id_BS,'missing_value','-9999.',/char

  ;;lonID_C = ncdf_vardef(id,'lonsC',[lnC_id],/float)
  ;;ncdf_attput,id,lonID_C,'units','degrees',/char
  ;;ncdf_attput,id,lonID_C,'long_name','longitudes in coarse grid',/char
  ;;latID_C = ncdf_vardef(id,'latsC',[ltC_id],/float)
  ;;ncdf_attput,id,latID_C,'units','degrees',/char
  ;;ncdf_attput,id,latID_C,'long_name','latitudes in coarse grid',/char

  ;;lonID_F = ncdf_vardef(id,'lonsF',[lnF_id],/double)
  ;;ncdf_attput,id,lonID_F,'units','degrees',/char
  ;;ncdf_attput,id,lonID_F,'long_name','longitudes in fine grid',/char
  ;;latID_F = ncdf_vardef(id,'latsF',[ltF_id],/double)
  ;;ncdf_attput,id,latID_F,'units','degrees',/char
  ;;ncdf_attput,id,latID_F,'long_name','latitudes in fine grid',/char

  ;;ncdf_attput,id,/global,'Title','Information of Broad Stratiform classified storms'
  ;;ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ;;ncdf_attput,id,/global,'rain_text1','rain calculated using ZR in Romatschke and Houze 2011'
  ;;ncdf_control,id,/endef                    ;; Put file in data mode:

  ;;ncdf_varput,id,freqF_id,freq_Full[*,*,3]                         ;;Write the data into file
  ;;ncdf_varput,id,freqC_id,freq_Core[*,*,3]

  ;;ncdf_varput,id,rainF_id_BS,rain_R11Full[*,*,3]
  ;;ncdf_varput,id,nR_F_id_BS,nRai_R11Full[*,*,3]
  ;;ncdf_varput,id,rain_id_BS,rain_R11Core[*,*,3]
  ;;ncdf_varput,id,nR_id_BS,nRai_R11Core[*,*,3]
  
  ;;ncdf_varput,id,lonID_C,lonsC
  ;;ncdf_varput,id,latID_C,latsC
  ;;ncdf_varput,id,lonID_F,lonsF
  ;;ncdf_varput,id,latID_F,latsF

  ;;ncdf_close,id
  ;;spawn,'zip -mjq '+name+'.zip '+name+'.nc'

  ;;;;save the accumulated rain matrices in Convective Cores 
  ;;name=path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/infoRainR11_Conv_EchoCores_'+type+'_v11s'
  ;;id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ;;ncdf_control,id,/fill                                  ;; Fill the file with default values
  ;;lnC_id = ncdf_dimdef(id,'lons_c',nlonsC)
  ;;ltC_id = ncdf_dimdef(id,'lats_c',nlatsC)
  ;;lnF_id = ncdf_dimdef(id,'lons_f',nlonsF)
  ;;ltF_id = ncdf_dimdef(id,'lats_f',nlatsF)
  ;;var_id = ncdf_dimdef(id,'echo_type',3)
  ;;freqF_id=ncdf_vardef(id,'freq_Full',[lnC_id,ltC_id,var_id],/long) ;; Define variables:
  ;;ncdf_attput,id,freqF_id,'long_name','Event pixel count for Entire storm',/char
  ;;ncdf_attput,id,freqF_id,'missing_value','-9999.',/char

  ;;freqC_id=ncdf_vardef(id,'freq_core',[lnC_id,ltC_id,var_id],/long) ;; Define variables:
  ;;ncdf_attput,id,freqC_id,'long_name','Event pixel count for cores',/char
  ;;ncdf_attput,id,freqC_id,'missing_value','-9999.',/char
  ;;;;** 
  ;;rainF_id_DC=ncdf_vardef(id,'rain_Full_DC',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rainF_id_DC,'long_name','Accumulated rain for Entire storm in DC',/char
  ;;ncdf_attput,id,rainF_id_DC,'missing_value','-9999.',/char

  ;;nR_F_id_DC=ncdf_vardef(id,'nRain_Full_DC',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_F_id_DC,'long_name','number of elements of rain for Entire storm in DC',/char
  ;;ncdf_attput,id,nR_F_id_DC,'missing_value','-9999.',/char

  ;;rain_id_DC=ncdf_vardef(id,'rain_Core_DC',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rain_id_DC,'long_name','Accumulated rain for DC cores',/char
  ;;ncdf_attput,id,rain_id_DC,'missing_value','-9999.',/char

  ;;nR_id_DC=ncdf_vardef(id,'nRain_Core_DC',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_id_DC,'long_name','number of elements of rain for DC cores',/char
  ;;ncdf_attput,id,nR_id_DC,'missing_value','-9999.',/char
  ;;;;** 
  ;;rainF_id_WC=ncdf_vardef(id,'rain_Full_WC',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rainF_id_WC,'long_name','Accumulated rain for Entire storm in WC',/char
  ;;ncdf_attput,id,rainF_id_WC,'missing_value','-9999.',/char

  ;;nR_F_id_WC=ncdf_vardef(id,'nRain_Full_WC',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_F_id_WC,'long_name','number of elements of rain for Entire storm in WC',/char
  ;;ncdf_attput,id,nR_F_id_WC,'missing_value','-9999.',/char

  ;;rain_id_WC=ncdf_vardef(id,'rain_Core_WC',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rain_id_WC,'long_name','Accumulated rain for WC cores',/char
  ;;ncdf_attput,id,rain_id_WC,'missing_value','-9999.',/char

  ;;nR_id_WC=ncdf_vardef(id,'nRain_Core_WC',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_id_WC,'long_name','number of elements of rain for WC cores',/char
  ;;ncdf_attput,id,nR_id_WC,'missing_value','-9999.',/char
  ;;;;** 
  ;;rainF_id_DW=ncdf_vardef(id,'rain_Full_DW',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rainF_id_DW,'long_name','Accumulated rain for Entire storm in DW',/char
  ;;ncdf_attput,id,rainF_id_DW,'missing_value','-9999.',/char

  ;;nR_F_id_DW=ncdf_vardef(id,'nRain_Full_DW',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_F_id_DW,'long_name','number of elements of rain for Entire storm in DW',/char
  ;;ncdf_attput,id,nR_F_id_DW,'missing_value','-9999.',/char

  ;;rain_id_DW=ncdf_vardef(id,'rain_Core_DW',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rain_id_DW,'long_name','Accumulated rain for DW cores',/char
  ;;ncdf_attput,id,rain_id_DW,'missing_value','-9999.',/char

  ;;nR_id_DW=ncdf_vardef(id,'nRain_Core_DW',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_id_DW,'long_name','number of elements of rain for DW cores',/char
  ;;ncdf_attput,id,nR_id_DW,'missing_value','-9999.',/char
  ;;;;**
  ;;lonID_C = ncdf_vardef(id,'lonsC',[lnC_id],/float)
  ;;ncdf_attput,id,lonID_C,'units','degrees',/char
  ;;ncdf_attput,id,lonID_C,'long_name','longitudes in coarse grid',/char
  ;;latID_C = ncdf_vardef(id,'latsC',[ltC_id],/float)
  ;;ncdf_attput,id,latID_C,'units','degrees',/char
  ;;ncdf_attput,id,latID_C,'long_name','latitudes in coarse grid',/char

  ;;lonID_F = ncdf_vardef(id,'lonsF',[lnF_id],/double)
  ;;ncdf_attput,id,lonID_F,'units','degrees',/char
  ;;ncdf_attput,id,lonID_F,'long_name','longitudes in fine grid',/char
  ;;latID_F = ncdf_vardef(id,'latsF',[ltF_id],/double)
  ;;ncdf_attput,id,latID_F,'units','degrees',/char
  ;;ncdf_attput,id,latID_F,'long_name','latitudes in fine grid',/char

  ;;ncdf_attput,id,/global,'Title','Information of Convective (DCC-WCC-DWCC) classified storms'
  ;;ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ;;ncdf_attput,id,/global,'rain_text1','rain calculated using ZR in Romatschke and Houze 2011'
  ;;ncdf_control,id,/endef                    ;; Put file in data mode:
  
  ;;ncdf_varput,id,freqF_id,freq_Full[*,*,[0,1,2]]                         ;;Write the data into file
  ;;ncdf_varput,id,freqC_id,freq_Core[*,*,[0,1,2]]

  ;;ncdf_varput,id,rainF_id_DC,rain_R11Full[*,*,0]
  ;;ncdf_varput,id,nR_F_id_DC,nRai_R11Full[*,*,0]
  ;;ncdf_varput,id,rain_id_DC,rain_R11Core[*,*,0]
  ;;ncdf_varput,id,nR_id_DC,nRai_R11Core[*,*,0]

  ;;ncdf_varput,id,rainF_id_WC,rain_R11Full[*,*,1]
  ;;ncdf_varput,id,nR_F_id_WC,nRai_R11Full[*,*,1]
  ;;ncdf_varput,id,rain_id_WC,rain_R11Core[*,*,1]
  ;;ncdf_varput,id,nR_id_WC,nRai_R11Core[*,*,1]

  ;;ncdf_varput,id,rainF_id_DW,rain_R11Full[*,*,2]
  ;;ncdf_varput,id,nR_F_id_DW,nRai_R11Full[*,*,2]
  ;;ncdf_varput,id,rain_id_DW,rain_R11Core[*,*,2]
  ;;ncdf_varput,id,nR_id_DW,nRai_R11Core[*,*,2]

  ;;ncdf_varput,id,lonID_C,lonsC
  ;;ncdf_varput,id,latID_C,latsC
  ;;ncdf_varput,id,lonID_F,lonsF
  ;;ncdf_varput,id,latID_F,latsF

  ;;ncdf_close,id
  ;;spawn,'zip -mjq '+name+'.zip '+name+'.nc'

  ;;save the Rain accumulation Matrices using Romatschke and Houze 2011!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;;name=path_out+'stats_class_v11s/RomHouze_2011/'+meses+'/infoRainR11_ShallowIsol_EchoCores_'+type+'_v11s'
  ;;id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ;;ncdf_control,id,/fill                                  ;; Fill the file with default values
  ;;lnC_id = ncdf_dimdef(id,'lons_c',nlonsC)
  ;;ltC_id = ncdf_dimdef(id,'lats_c',nlatsC)
  ;;lnF_id = ncdf_dimdef(id,'lons_f',nlonsF)
  ;;ltF_id = ncdf_dimdef(id,'lats_f',nlatsF)

  ;;freqF_id=ncdf_vardef(id,'freq_Full',[lnC_id,ltC_id],/long) ;; Define variables:
  ;;ncdf_attput,id,freqF_id,'long_name','pixel count for full storm',/char
  ;;ncdf_attput,id,freqF_id,'missing_value','-9999.',/char

  ;;freqC_id=ncdf_vardef(id,'freq_core',[lnC_id,ltC_id],/long) ;; Define variables:
  ;;ncdf_attput,id,freqC_id,'long_name','pixel count for core storm',/char
  ;;ncdf_attput,id,freqC_id,'missing_value','-9999.',/char
  ;;;;**    
  ;;rainF_id_SH=ncdf_vardef(id,'rain_Full_SH',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rainF_id_SH,'long_name','Accumulated rain for full storm in SH',/char
  ;;ncdf_attput,id,rainF_id_SH,'missing_value','-9999.',/char

  ;;nR_F_id_SH=ncdf_vardef(id,'nRain_Full_SH',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_F_id_SH,'long_name','number of elements of rain for full storm in SH',/char
  ;;ncdf_attput,id,nR_F_id_SH,'missing_value','-9999.',/char

  ;;rain_id_SH=ncdf_vardef(id,'rain_Core_SH',[lnF_id,ltF_id],/float) ;; Define variables:
  ;;ncdf_attput,id,rain_id_SH,'long_name','Accumulated rain for SH cores',/char
  ;;ncdf_attput,id,rain_id_SH,'missing_value','-9999.',/char

  ;;nR_id_SH=ncdf_vardef(id,'nRain_Core_SH',[lnF_id,ltF_id],/long) ;; Define variables:
  ;;ncdf_attput,id,nR_id_SH,'long_name','number of elements of rain for SH cores',/char
  ;;ncdf_attput,id,nR_id_SH,'missing_value','-9999.',/char

  ;;lonID_C = ncdf_vardef(id,'lonsC',[lnC_id],/float)
  ;;ncdf_attput,id,lonID_C,'units','degrees',/char
  ;;ncdf_attput,id,lonID_C,'long_name','longitudes in coarse grid',/char
  ;;latID_C = ncdf_vardef(id,'latsC',[ltC_id],/float)
  ;;ncdf_attput,id,latID_C,'units','degrees',/char
  ;;ncdf_attput,id,latID_C,'long_name','latitudes in coarse grid',/char

  ;;lonID_F = ncdf_vardef(id,'lonsF',[lnF_id],/double)
  ;;ncdf_attput,id,lonID_F,'units','degrees',/char
  ;;ncdf_attput,id,lonID_F,'long_name','longitudes in fine grid',/char
  ;;latID_F = ncdf_vardef(id,'latsF',[ltF_id],/double)
  ;;ncdf_attput,id,latID_F,'units','degrees',/char
  ;;ncdf_attput,id,latID_F,'long_name','latitudes in fine grid',/char

  ;;ncdf_attput,id,/global,'Title','Information of Shallow Isolated storms'
  ;;ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ;;ncdf_attput,id,/global,'rain_text1','rain calculated using ZR in Romatschke and Houze 2011'
  ;;ncdf_control,id,/endef                    ;; Put file in data mode:

  ;;ncdf_varput,id,freqF_id,freq_Full[*,*,4]                         ;;Write the data into file
  ;;ncdf_varput,id,freqC_id,freq_Core[*,*,4]

  ;;ncdf_varput,id,rainF_id_SH,rain_R11Full[*,*,4]
  ;;ncdf_varput,id,nR_F_id_SH,nRai_R11Full[*,*,4]
  ;;ncdf_varput,id,rain_id_SH,rain_R11Core[*,*,4]
  ;;ncdf_varput,id,nR_id_SH,nRai_R11Core[*,*,4]

  ;;ncdf_varput,id,lonID_C,lonsC
  ;;ncdf_varput,id,latID_C,latsC
  ;;ncdf_varput,id,lonID_F,lonsF
  ;;ncdf_varput,id,latID_F,latsF

  ;;ncdf_close,id
  ;;spawn,'zip -mjq '+name+'.zip '+name+'.nc'
  ;;*********************************************************************************************************************

  ;;*********************************************************************************************************************
  ;;Save the cfads for BSR storms - check for directories
  if file_test(path_out+'stats_class_v11s/Cfad',/directory) eq 0 then file_mkdir,path_out+'stats_class_v11s/Cfad'
  if file_test(path_out+'stats_class_v11s/Cfad/'+meses,/directory) eq 0 then  $
     file_mkdir,path_out+'stats_class_v11s/Cfad/'+meses
  
  name=path_out+'stats_class_v11s/Cfad/'+meses+'/infoCfad_Stra_EchoCores_'+type+'_v11s'
  id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ncdf_control,id,/fill                                  ;; Fill the file with default values
  alt_id = ncdf_dimdef(id,'alts_CFAD',nlevels)
  ref_id = ncdf_dimdef(id,'refl_CFAD',n_refls)

  cfadF_id=ncdf_vardef(id,'CFAD_Full',[ref_id,alt_id],/long)
  ncdf_attput,id,cfadF_id,'long_name','CFAD count for Entire storm',/char
  ncdf_attput,id,cfadF_id,'missing_value','-9999.',/char

  cfadC_id=ncdf_vardef(id,'CFAD_Core',[ref_id,alt_id],/long)
  ncdf_attput,id,cfadC_id,'long_name','CFAD count for cores',/char
  ncdf_attput,id,cfadC_id,'missing_value','-9999.',/char

  ncdf_attput,id,/global,'Title','CFADs for BSR storms'
  ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ncdf_control,id,/endef                    ;; Put file in data mode:

  ncdf_varput,id,cfadF_id,CFAD_Full[*,*,3]
  ncdf_varput,id,cfadC_id,CFAD_Core[*,*,3]

  ncdf_close,id
  spawn,'zip -mjq '+name+'.zip '+name+'.nc'

  ;;Save the cfads for Convective storms
  name=path_out+'stats_class_v11s/Cfad/'+meses+'/infoCfad_Conv_EchoCores_'+type+'_v11s'
  id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ncdf_control,id,/fill                                  ;; Fill the file with default values
  alt_id = ncdf_dimdef(id,'alts_CFAD',nlevels)
  ref_id = ncdf_dimdef(id,'refl_CFAD',n_refls)
  var_id = ncdf_dimdef(id,'echo_type',3)
  
  cfadF_id=ncdf_vardef(id,'CFAD_Full',[ref_id,alt_id,var_id],/long)
  ncdf_attput,id,cfadF_id,'long_name','CFAD count for Entire storm',/char
  ncdf_attput,id,cfadF_id,'missing_value','-9999.',/char

  cfadC_id=ncdf_vardef(id,'CFAD_Core',[ref_id,alt_id,var_id],/long)
  ncdf_attput,id,cfadC_id,'long_name','CFAD count for cores',/char
  ncdf_attput,id,cfadC_id,'missing_value','-9999.',/char

  ncdf_attput,id,/global,'Title','CFADs for Convective storms'
  ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ncdf_control,id,/endef                    ;; Put file in data mode:

  ncdf_varput,id,cfadF_id,CFAD_Full[*,*,[0,1,2]]
  ncdf_varput,id,cfadC_id,CFAD_Core[*,*,[0,1,2]]

  ncdf_close,id
  spawn,'zip -mjq '+name+'.zip '+name+'.nc'

  ;;Save the CFAD accumulation matrices for Shallow Isolated storms
  name=path_out+'stats_class_v11s/Cfad/'+meses+'/infoCfad_ShallowIsol_EchoCores_'+type+'_v11s'
  id = ncdf_create(name+'.nc',/clobber)                  ;;Create NetCDF file
  ncdf_control,id,/fill                                  ;; Fill the file with default values
  alt_id = ncdf_dimdef(id,'alts_CFAD',nlevels)
  ref_id = ncdf_dimdef(id,'refl_CFAD',n_refls)

  cfadF_id=ncdf_vardef(id,'CFAD_Full',[ref_id,alt_id],/long)
  ncdf_attput,id,cfadF_id,'long_name','CFAD count for Entire storm',/char
  ncdf_attput,id,cfadF_id,'missing_value','-9999.',/char

  cfadC_id=ncdf_vardef(id,'CFAD_Core',[ref_id,alt_id],/long)
  ncdf_attput,id,cfadC_id,'long_name','CFAD count for cores',/char
  ncdf_attput,id,cfadC_id,'missing_value','-9999.',/char

  ncdf_attput,id,/global,'Title','CFADs for Shallow Isolated storms'
  ncdf_attput,id,/global,'source','created using IDL by S.Brodzik UW, Seattle'
  ncdf_control,id,/endef                    ;; Put file in data mode:
  
  ncdf_varput,id,cfadF_id,CFAD_Full[*,*,4]
  ncdf_varput,id,cfadC_id,CFAD_Core[*,*,4]

  ncdf_close,id
  spawn,'zip -mjq '+name+'.zip '+name+'.nc'
  ;;*********************************************************************************************************************

  print,'start='+comienzo
  print,'end='+systime(/utc)

end

function deg2km,deg,lon,lat
  ;;**Function to compute the distance in km from a distance in degrees for a region of
  ;;**the Earth centered at the lat-lon coordinates
  ;;kms[0] - Zonal distance (Equivalent of 1°) in kilometers at given Lon, Lat
  ;;kms[1] - meridional distance (Equivalent of 1° )in Kiometers at given Lon, Lat and
  ;;Re=6378.137   ;;this is the earth's radius at the equator   
  Re=6371.      ;;mean earth radius (this is equivalent to the matlab function)
  kms=fltarr(2)
  kms[0]=Re/(180.0/!pi)
  kms[1]=(Re/(180.0/!pi))*cos(lat*(!pi/180))
  kms=deg*kms
  return,kms
end

PRO undefine, varname  
  tempvar = SIZE(TEMPORARY(varname))
END


