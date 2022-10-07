from math import *

def dist_and_bearing_between_two_latlons(lat1,lon1,lat2,lon2):

    #lat1 = 0
    #lat2 = -1
    #lon1 = 30
    #lon2 = 30

    #lat1 = lat_dow
    #lon1 = lon_dow
    #lat2 = lat_value
    #lon2 = lon_value
    #print(lat1,lon1,lat2,lon2)
    
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    Base = 6371 * c

    Bearing = atan2(sin(lon2-lon1)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon2-lon1))
    Bearing = degrees(Bearing)
    Bearing = (Bearing + 360) % 360

    #dist = Base
    #bearing = Bearing
    #print(Base)
    #print(Bearing)

    return Base, Bearing
