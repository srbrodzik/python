#!/bin/csh

if ($#argv != 1) then
  echo Usage: $0 region
  exit 1
endif

set region = $1
echo $region

./convert_info_rain_nsr_bsr_grd_cdf_to_cdf.py $region str
./convert_info_rain_nsr_bsr_grd_cdf_to_cdf.py $region mod
./convert_info_rain_nsr_conv_grd_cdf_to_cdf.py $region str
./convert_info_rain_nsr_conv_grd_cdf_to_cdf.py $region mod
./convert_info_rain_nsr_shi_grd_cdf_to_cdf.py $region
./convert_monthly_class_tab_ascii_to_cdf.py $region str
./convert_monthly_class_tab_ascii_to_cdf.py $region mod
./convert_cfad_dbz_grd_cdf_to_cdf.py $region str
./convert_cfad_dbz_grd_cdf_to_cdf.py $region mod
./convert_cfad_dm_grd_cdf_to_cdf.py $region str
./convert_cfad_dm_grd_cdf_to_cdf.py $region mod
