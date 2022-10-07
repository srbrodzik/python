# Copied this script from pyart-users post and thought it might be useful down
# the road.  Looks like it's reading nexrad data and interpolating it.  SRB

## START GRIDDING THE DATA
# read in the data
outputdir = '/home/kkossen/radar_files/' + folder_name.replace('_gz', '') +
'_gridded_filtered/'
filename = 'K*_*_V06*'
path = datadir + filename
files = sorted(glob.glob(path))
for file in files:
radar = pyart.io.read(file)
print("")
 
#make a string that is the output file name
outfilename = file.replace(datadir, outputdir)
 
# mask out last 10 gates of each ray, this removes the "ring" around the
radar.
radar.fields['reflectivity']['data'][:, -10:] = np.ma.masked
 
# exclude masked gates from the gridding
gatefilter = pyart.filters.GateFilter(radar)
gatefilter.exclude_transition()
gatefilter.exclude_masked('reflectivity')
 
print("Creating grid...")
 
# perform Cartesian mapping, limit to the reflectivity field.
grid2 = pyart.map.grid_from_radars(
(radar,), gatefilters=(gatefilter, ),
grid_shape=(50, 500, 500),
grid_limits=((0, 15000), (-250000.0, 250000.0), (-250000.0,
350000.0)),
fields=['reflectivity']
#roi_func = 'dist',
#z_factor = 0.1,
#xy_factor = 0.02
)
print("Grid created")
 
#write netCDF for grid2
pyart.io.write_grid('%s_filtered_grid', grid2, format = 'NETCDF4')
print("NetCDF file created")
 
#recall made netCDF files
path = outputdir + 'K*_*_filtered_grid'
grids = sorted(glob.glob(path))
 
for file in grids:
grid = pyart.io.read_grid(file)
# Do whatever you want with grid (numpy array)
 
 
