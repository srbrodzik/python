import math

def compute_latlon(lat_in,lon_in,brng_in,d):

    R = 6378.1 #Radius of the Earth
    #brng = math.radians(90.0) #Bearing 
    #d = 15 #Distance in km

    #TEST CASE
    #lat1  52.20472
    #lon1  0.14056
    #brng  90.0 deg
    #d     15
    #lat2  52.20444 - the lat result I'm hoping for
    #lon2  0.36056 - the long result I'm hoping for.

    # convert degrees to radians
    lat1 = math.radians(lat_in)
    lon1 = math.radians(lon_in)
    brng = math.radians(brng_in)

    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
                      math.cos(lat1)*math.sin(d/R)*math.cos(brng))

    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
                             math.cos(d/R)-math.sin(lat1)*math.sin(lat2))
    
    lat_out = math.degrees(lat2)
    lon_out = math.degrees(lon2)

    return lat_out,lon_out
    
    #print(lat2)
    #print(lon2)
