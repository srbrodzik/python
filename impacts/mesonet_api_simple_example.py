# Copied from
# https://github.com/joejoezz/blog/blob/master/2018/mesonet_api_simple_example.py
# which is based on
# https://synopticlabs.org/synoptic (brodzik@uw.edu; xxxxxsynop)
#   Default public token: 8c150f37a5ba47ad929c0a24180e877c
#   Private API key: Ky46jvXMfpPhKKqWBpGwRa15e1y19MfAapEK8ESeHZ

import urllib
import json
import pandas as pd
import matplotlib.pyplot as plt

# Specify request parameters (as strings)
token = # use your own token

# Create API query string
args = {
    'start':'201806010000',
    'end':'201806070000',
    'obtimezone':'UTC',
    'vars':'air_temp',
    'stids':'KSEA',
    'units':'temp|F',
    'token':token
}

apiString = urllib.urlencode(args)
url = "http://api.mesowest.net/v2/stations/timeseries"
fullUrl = '{}?{}'.format(url,apiString)

# Open the URL and convert the returned JSON into a dictionary
response = urllib.urlopen(fullUrl)
responseDict = json.loads(response.read())

# Isolate the time and temperature from the response dictionary
dateTime = responseDict['STATION'][0]['OBSERVATIONS']['date_time']
airT = responseDict['STATION'][0]['OBSERVATIONS']['air_temp_set_1']
ksea = pd.Series(airT,index=pd.to_datetime(dateTime))

# Retain only the hourly observations
ksea = ksea.where(ksea.index.minute == 53).dropna()

# Plotting code
fig,ax = plt.subplots()
ksea.plot(ax=ax)
plt.title('KSEA Hourly Temperature ($^{\circ}$F)')

# Clean up the plot a bit
from matplotlib.dates import DayLocator, DateFormatter
ax.xaxis.set_major_locator(DayLocator())
ax.xaxis.set_major_formatter(DateFormatter('%d-%b'))
ax.set_ylabel('Temperature ($^{\circ}$F)')
ax.grid()

plt.savefig('temp2.png')
