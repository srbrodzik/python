import math

# For Gulf
midlat = (32.5-25.)/2 + 25.
midlon = -100. + (-75. + 100.)/2.

lat1 = midlat
lat2 = midlat
lon1 = -100.
lon2 = -75.

lat1 = 25
lat2 = 32.5
lon1 = midlon
lon2 = midlon


# For Plains
midlat = (55. - 35.)/2 + 35.
midlon = -105. + (-90. + 105.)/2.

lat1 = midlat
lat2 = midlat
lon1 = -105.
lon2 = -90.

lat1 = 35
lat2 = 55
lon1 = midlon
lon2 = midlon

def convert_latLonDist_to_km(lat1, lon1, lat2, lon2):
    R = 6378.137; # Radius of earth in KM
    dLat = lat2 * math.pi / 180. - lat1 * math.pi / 180.
    dLon = lon2 * math.pi / 180. - lon1 * math.pi / 180.
    a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(lat1 * math.pi / 180.) * math.cos(lat2 * math.pi / 180.) * math.sin(dLon/2) * math.sin(dLon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = R * c
    return d # kilometers
