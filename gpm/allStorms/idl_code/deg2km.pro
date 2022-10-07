function deg2km,deg,lon,lat
  ;;**Function to compute the distance in km from a distance in degrees for a region of
  ;;**the Earth centered at the lat-lon coordinates
  ;;kms[0] - Zonal distance (Equivalent of 1°) in kilometers at given Lon, Lat
  ;;kms[1] - meridional distance (Equivalent of 1° )in Kiometers at given Lon, Lat and
  ;;Re=6378.137   ;;this is the earth's radius at the equator   
  Re=6371.      ;;mean earth radius (this is equivalent to the matlab function)
  kms=fltarr(2)
  kms[0]=Re/(180.0/!pi)
  kms[1]=(Re/(180.0/!pi))*cos(lat*(!pi/180))
  kms=deg*kms
  return,kms
end
