pro run_allStorms_wmp_v11s_v06

  ;;  .r /home/disk/shear2/brodzik/IDL/gpm/allStorms/LandRuns/run_allStorms_wmp_v11s_v06

  ;;resolve_routine,'findStorm'   ;;.r /home/disk/meso-home/mzuluaga/MesoAmerica/classify/id_storms/findStorm
  resolve_routine,'findStorm'   ;;  .r /home/disk/shear2/brodzik/IDL/gpm/findStorm
  resolve_routine,'findStormNew';;  .r /home/disk/shear2/brodzik/IDL/gpm/findStormNew
  
  resolve_routine,'allStorms_wmp_v11s_v06'   ;;  .r /home/disk/shear2/brodzik/IDL/gpm/allStorms/LandRuns/allStorms_wmp_v11s_v06

  ;;years=['2014','2015','2016','2017','2018','2019'] 
  years=['2016'] 
  ;;month=['01','02','03','04','05','06','07','08','09','10','11','12']  
  month=['09']
  
  region='WMP'
  
  ;;print,years,month

  for yy=0,n_elements(years)-1 do begin

     ;;if years[yy] eq '2018' then month = ['12'] else  $
        ;;if years[yy] eq '2019' then month = ['01','02'] else $
           ;;month=['01','02','03','04','05','06','07','08','09','10','11','12']  

     for mm=0,n_elements(month)-1 do begin
        typeFile=month[mm]+'_'+years[yy]+'_'+region
        allStorms_wmp_v11s_v06,type=typeFile
     endfor

  endfor 
   
  stop
end


;; directory
;;path_in ='/home/disk/bob/gpm/wmp_ku/classify/ex_data_v06/'
;;path_out='/home/disk/bob/gpm/wmp_ku/classify/class_data_v06/'

