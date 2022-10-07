
# coding: utf-8

# In[1]:

# Block 0: Documentation

print('Program to list and download ABI data files from Amazon Web Services (AWS), AMS Short Course, March 18, 2021\n')
print('Version 1.0, February 8, 2021\n')
print('Written by Dr. Amy Huff (IMSG at NOAA/NESDIS/STAR) and Ryan Theurer (GVT LLC at NOAA/NESDIS/STAR)\n')
print('For questions contact Dr. Huff: amy.huff@noaa.gov\n')
print('This program accesses the ABI data archive on AWS and lists/downloads data files for a given date and time period specified by the user.\n')
print('Block 1 imports modules and libraries, and blocks 2-3 are functions that require no input from the user; there is no visible output from these blocks. In block 4, the user enters settings and obtains output.')


# In[21]:

# Block 1: Import modules and libraries

# Library to perform array operations
import numpy as np

# Module to interface with s3 (AWS)
import s3fs

# Module for manipulating dates and times
import datetime

# Module to access files in the directory
import os


# In[34]:

# Block 2: Find Julian day from given year/month/day
# "year", "month", and "day" are global variables set in final block

def julian(year, month, day):
    calendar = datetime.datetime(year, month, day)
    julian_day = calendar.strftime('%j')
    
    return julian_day


# In[35]:

# Block 3: Create array containing ABI data file names for given satellite/product and date/time period
# "year", "month", "day", "start", "end", "satellite", and "product" are global variables set in final block

def aws_list(year, month, day, start, end, satellite, product):
  
    # Access AWS using anonymous credentials
    aws = s3fs.S3FileSystem(anon=True)

    # Make a list of all data files encompassing given date and start/end hours
    julian_day = julian(int(year), int(month), int(day))
    start_time = start[0:2]
    end_time = end[0:2]
    hour_range = range(int(start_time), int(end_time) + 1)
    final_list = []
    for i in hour_range:
        ii = str(i)
        if len(ii) < 2:
            ii = '0'+ii
        hour_files = aws.ls('noaa-goes' + str(satellite) + '/' + product + '/' + year + '/' + julian_day + '/' + ii + '/')
        final_list.extend(hour_files) 
        all_hours = np.array(final_list)

    # Extract list of data files for specified period set by start/end times
    data = []
    # List file names
    for i in all_hours:
        if i[-42:-38] >= start and i[-42:-38] <= end:
            data.append(i)
        else:
            continue

    return data


# In[63]:

# Block 4: Enter user settings and list ABI data files, with option to save files locally

# Data file saving settings
save = 'yes'  # Option to save data files: 'yes' (save to "file_path" directory) or 'no' (list file names only)
#save_path = os.getcwd() + '/data/'  # Directory where data files will be saved
save_path = '/home/disk/monsoon/precip/raw/satellite/goes17/FullDisk/Channel09/'  # Directory where data files will be saved

# Satellite and product settings
satellite = 17  # GOES-East = 16, GOES-West = 17
#product = 'ABI-L2-AODC'  # ABI product name abbreviation; see list at https://docs.opendata.aws/noaa-goes16/cics-readme.html
product = 'ABI-L1b-RadF'  # ABI product name abbreviation; see list at https://docs.opendata.aws/noaa-goes16/cics-readme.html
channel = 'C09'

# Day and time settings
year = '2021'    # 4-digit year (e.g., 2021)
#month = '7'    # 1- or 2-digit month (e.g., Feb = 2, Oct = 10)
month = '8'    # 1- or 2-digit month (e.g., Feb = 2, Oct = 10)
#days = ['15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']    # 1- or 2- digit day (e.g., 1, 25)
days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19']    # 1- or 2- digit day (e.g., 1, 25)
#start = '1200'    # 4-digit observation start time in UTC, no colon (e.g. 20:00 UTC = 2000)
start = '0000'    # 4-digit observation start time in UTC, no colon (e.g. 20:00 UTC = 2000)
end = '2359'    # 4-digit observation end time in UTC, no colon (e.g. 20:35 UTC = 2035)
#start = '0000'    # 4-digit observation start time in UTC, no colon (e.g. 20:00 UTC = 2000)
#end = '0300'    # 4-digit observation end time in UTC, no colon (e.g. 20:35 UTC = 2035)

#########################################################################################################################

if __name__ == "__main__":

    for day in days:
    
        # Query AWS and list filenames matching entered settings
        data = aws_list(year, month, day, start, end, satellite, product)
        if len(data) > 0:
            for i in data:
                if channel in i:
                    print(i.split('/')[-1])
        else:
            print('No files retrieved.  Check settings and try again.')

        # Downlad and save data files to specfied directory
        if save == 'yes':
            if not os.path.isdir(save_path):
                os.makedirs(save_path)
            aws = s3fs.S3FileSystem(anon=True)
            for i in data:
                if channel in i:
                    print(i)
                    aws.get(i, save_path + i.split('/')[-1])
            print('Download complete!')
        else:
            pass


# In[ ]:



