#!/bin/bash

# Script to run one day of GPM orbits
# provided by Dr. Rachael Kroodsma as example with hardcoded locations.
# then modified to accept input parameters that provid the year, month,
# days and the output directory in which products are made.

# The python program sends in the appropriate year, day and month as well 
# as the output directory. It
# executes 11 days in parallel as that is the max I can do with our matlab.

YYYY="$1"
MM="$2"
Days="$3"
#Days=`echo {01..10}`

dirOUT=$4
CCdir=$HOME/stacyKU

#################################################################
# Generally the Days variable will be a single day of the month #
# to process. #
#################################################################
for DD in $Days
do
 dirKu=/gpmdata/$YYYY/$MM/$DD/radar/
 $CCdir/getGPM_Ku_interp_SAM_v3_compile $dirOUT $dirKu
 $CCdir/getGPM_Ku_interp_ASIA_v3_compile $dirOUT $dirKu
done

   
