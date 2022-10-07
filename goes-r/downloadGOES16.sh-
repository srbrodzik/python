#!/bin/bash
 
#########################################
# LICENSE
#Copyright (C) 2012 Dr. Marcial Garbanzo Salas

#This program is free software: you can redistribute it and/or modify it under the terms of the
#GNU General Public License as published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
#even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program. If not,
#see http://www.gnu.org/licenses/.
#########################################
 
#########################################
# AUTHOR
# This program was created at the University of Costa Rica (UCR)
# It is intended as a tool for meteorology students to obtain data from GOES16
# but it can be used by operational and research meteorology.
#########################################
 
#########################################
# Warning: This program can download a LARGE amount of information
# and this can cause problems with limited bandwidth networks or
# computers with low storage capabilities.
#########################################
 
#########################################
# CLEANING FROM PREVIOUS RUNS
#
rm DesiredData.txt
rm FinalData.txt
rm FullList.txt
#########################################
 
echo "GOES16 ABI data downloader"
 
#########################################
# CONFIGURATION
#
# YEAR OF INTEREST
YEARS='2021'
 
# DAYS OF THE YEAR
# Can use this link to find out: https://www.esrl.noaa.gov/gmd/grad/neubrew/Calendar.jsp
# Example: 275 for October 2nd, 2017
# NOTE: There is only about 60 days previous to the current date available
DAYS="058"

#HOURS='00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
HOURS='09 10 11 12 13 14'
 
# CHANNELS (micrometers)
# Example: CHANNELS='C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14 C15 C16'
# C01     0.47     Visible - Blue
# C02     0.64     Visible - Red
# C03     0.86     Near IR - Vegetation
# C04     1.37     Near IR - Cirrus
# C05     1.6      Near IR - Snow/Ice
# C06     2.2      Near IR - Cloud Ice
# C07     3.9      Near IR - Shortwave
# C08     6.2      IR - Upper-Level WV
# C09     6.9      IR - Mid-Level WV
# C10     7.3      IR - Lower-Level WV
# C11     8.4      IR - Cloud Top Phase
# C12     9.6      IR - Ozone
# C13     10.3     IR - Clean
# C14     11.2     IR - Standard
# C15     12.3     IR - Dirty
# C16     13.3     IR - CO2
CHANNELS='C01 C02 C03 C08 C09 C10 C13'
#CHANNELS='C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14 C15 C16'

# ABI PRODUCTS
# For a description look into:
# https://aws.amazon.com/public-datasets/goes/
# and
# http://edc.occ-data.org/goes16/getdata/
# Example: PRODUCTS='L1b-RadC L1b-RadF L1b-RadM L2-CMIPC L2-CMIPF L2-CMIPM L2-MCMIPC L2-MCMIPF L2-MCMIPM'
# ABI = Advanced Baseline Imager
# Radiance data (Level 1b)
#     L1b-RadC = CONUS
#     L1b-RadF = Full Disk
#     L1b-RadM = Mesoscale
# CMI = Cloud and Moisture Imagery (Level 2)
#     L2-CMIPC = CONUS
#     L2-CMIPF = Full Disk
#     L2-CMIPM = Mesoscale
# MCM = Multi-Band Cloud & Moisture Imagery (Level 2)
#     L2-MCMIPC = CONUS
#     L2-MCMIPF = Full Disk
#     L2-MCMIPF = Mesoscale
PRODUCTS='L1b-RadC'
#PRODUCTS='L2-CMIPF'
#PRODUCTS='L1b-RadM'

# SECTOR(S)
# If looking for mesoscale sector data, indicate with sector you want (M1 or M2)
# Example: SECTORS='M1 M2'
SECTORS=''
#########################################
 
#########################################
# Get list of remote files available
# PART 1. Obtain full list of files
#
for PRODUCT in $PRODUCTS; do
    for YEAR in $YEARS; do
	for DAY in $DAYS; do
	    for HOUR in $HOURS; do

		aws s3 --no-sign-request ls --recursive noaa-goes16/ABI-$PRODUCT/$YEAR/$DAY/$HOUR/ | awk '{print $3";"$4}' >> FullList.txt

	    done
	done
    done
done
 
#
# PART 2. Select only desired channels and sectors
for CHANNEL in $CHANNELS; do
    grep $CHANNEL FullList.txt >> DesiredData.txt
done

if [ $PRODUCTS == 'L1b-RadM' ]; then
    for SECTOR in $SECTORS; do
	grep $SECTOR DesiredData.txt >> FinalData.txt
    done
fi
#########################################
 
#########################################
# DOWNLOAD
#

# NOTE: If looking for meso sector data, use FinalData.txt in next line
#       of code instead of DesiredData.txt
for x in $(cat DesiredData.txt);
do
    SIZE=$(echo $x | cut -d";" -f1)
    FULLNAME=$(echo $x | cut -d";" -f2)
    NAME=$(echo $x | cut -d"/" -f5)
 
    echo "Processing file $NAME of size $SIZE"
    if [ -f $NAME ]; then
	echo "This file exists locally"
	LOCALSIZE=$(du -sb $NAME | awk '{ print $1 }')
	if [ $LOCALSIZE -ne $SIZE ]; then
	    echo "The size of the file is not the same as the remote file. Downloading again..."
	    aws s3 --no-sign-request cp s3://noaa-goes16/$FULLNAME ./
	else
	    echo "The size of the file matches the remote file. Not downloading it again."
	fi
    else
	echo "This file does not exists locally, downloading..."
	aws s3 --no-sign-request cp s3://noaa-goes16/$FULLNAME ./
    fi
 
done
#########################################
 
echo Program ending.
