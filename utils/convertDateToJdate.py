# Input: String YYYYMMDD
# Output: String YYYYjjj

import datetime
import sys

if len(sys.argv) != 2:
    print('Usage: {} [YYYYMMDD]'.format(sys.argv[0]))
    sys.exit()
else:
    date = sys.argv[1]

# calculate julian day
dt = datetime.datetime.strptime(date,'%Y%m%d')
yearStr = dt.strftime('%Y')
tt = dt.timetuple()
jday = tt.tm_yday
jdate = yearStr+str(jday).zfill(3)

return jdate
