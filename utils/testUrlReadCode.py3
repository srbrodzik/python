#!/usr/bin/python3

from getFilesFromProtUrl import get_flist

url = 'https://doppler.somas.stonybrook.edu/IMPACTS/BNL/mrrpro2white/images/quicklooks'
user = 'DataAccess'
pwd = 'WinterAtSBU'
suffix = 'png'

flist = get_flist(url,user,pwd,suffix)
print(flist)

