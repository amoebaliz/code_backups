# load the modules
from pydap.client import open_url
import pprint
import matplotlib.pylab as plt
import numpy as np
import datetime as dt
import pyroms
import find_geographical_pt as findpt
import netCDF4 as nc

my_target_grid='CORAL'
#my_target_grid='NWA'
#my_target_grid='CCS1'
I_dont_know_my_subset=False
my_workspace = '/t3/workdir/raphael/hycom_4_coral/'

# create an url for the HYCOM dataset
urlroot='http://tds.hycom.org/thredds/dodsC/GLBa0.08/'
myurl = urlroot + 'expt_90.6'
# import the dataset
dataset = open_url(myurl)

# optional : print the variables
#pprint.pprint( dataset.keys() )

#-------------------------------------------------------------------------------
# finding the subset we need for our particular target ROMS grid
# this is a bit long so it is better to run it once then copy the result into
# the else 
if I_dont_know_my_subset is True:
	# get lon/lat for ROMS target grid
	trg_grd = pyroms.grid.get_ROMS_grid(my_target_grid)
	trg_lon_min = trg_grd.hgrid.lon_vert.min()
	trg_lon_max = trg_grd.hgrid.lon_vert.max()
	trg_lat_min = trg_grd.hgrid.lat_vert.min()
	trg_lat_max = trg_grd.hgrid.lat_vert.max()

	print 'Target grid ' + my_target_grid + ' longitude span from ', trg_lon_min , '->', trg_lon_max
	print 'Target grid ' + my_target_grid + ' latitude span from ' , trg_lat_min , '->', trg_lat_max

	# get global lon/lat from hycom
	lon_hycom = dataset['Longitude']['Longitude'][:]
	lat_hycom = dataset['Latitude']['Latitude'][:]

	# find what are the range of indices in the global that we need in our area
	trg = findpt.find_geo(lon_hycom,lat_hycom)
	jmin, imin = trg(trg_lon_min,trg_lat_min)
	jmax, imax = trg(trg_lon_max,trg_lat_max)

	# because we're paranoid, let's take N more points on each side
	npts=12
	imin=imin-npts ; jmin=jmin-npts ; imax=imax+npts ; jmax=jmax+npts

	print 'Write down this so that you dont have to run this step everytime'
	print 'Lower Left corner (imin,jmin) = ', imin,jmin
	print 'Upper Right corner (imax,jmax) = ', imax,jmax

	ssh_test = dataset['ssh']['ssh'][0,jmin:jmax,imin:imax].squeeze()
	lon_test = dataset['Longitude']['Longitude'][jmin:jmax,imin:imax].squeeze()
	lat_test = dataset['Latitude']['Latitude'][jmin:jmax,imin:imax].squeeze()
	sshplt = np.ma.masked_values(ssh_test,ssh_test.max())
	# verification plot
	print 'Plotting for verification...'
	plt.figure() ; plt.contourf(lon_test,lat_test,sshplt,35) ; plt.colorbar()  
	plt.plot(trg_grd.hgrid.lon_vert[0,:],trg_grd.hgrid.lat_vert[0,:],'ko')
	plt.plot(trg_grd.hgrid.lon_vert[-1,:],trg_grd.hgrid.lat_vert[-1,:],'ko')
	plt.plot(trg_grd.hgrid.lon_vert[:,0],trg_grd.hgrid.lat_vert[:,0],'ko')
	plt.plot(trg_grd.hgrid.lon_vert[:,-1],trg_grd.hgrid.lat_vert[:,-1],'ko')
	plt.title('Sample ssh for grid ' + my_target_grid)
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.show()
	exit()
else:
	# I know what my subset is
	# CORAL
	if my_target_grid == 'CORAL':
		imin=245 ; jmin=1196 ; imax=1188 ; jmax=1870
	else:
		'grid indices not know, run me with I_dont_know_my_subset=True'


#-------------------------------------------------------------------------------
# Now that I know what to extract, let's proceed

# time vector in hycom does not have the same reference
# so we correct this first
time_vector = dataset['MT'].data[:]
ref_string  = dataset['MT'].attributes['units']
fmt = 'days since %Y-%m-%d %H:%M:%S'
ref_hycom = dt.datetime.strptime(ref_string,fmt)
# the person that wrote the unit attribute most likely
# made a mistake thinking that 1900-12-31 00:00:00 was 
# the end of year but it is the day before
ref_hycom = ref_hycom + dt.timedelta(days=1)
# roms uses jan 1 1900 as a reference usually
ref_roms = dt.datetime(1900,1,1,0,0)
# correcting the time vector
diff = (ref_hycom - ref_roms).days
time_vector_corrected = time_vector + diff

# getting invariant data (lon/lat/depth) for our region
lon_out   = dataset['Longitude']['Longitude'][jmin:jmax,imin:imax]
lat_out   = dataset['Latitude']['Latitude'][jmin:jmax,imin:imax]
depth_out = dataset['Depth'].data[:]

# size of arrays
ny = lat_out.shape[0]
nx = lat_out.shape[1]
nz = depth_out.shape[0]

# number of frames to process
nf = len(time_vector_corrected)

for kf in np.arange(nf):
	# retrieve data from opendap
	zeta = dataset['ssh']['ssh'][kf,jmin:jmax,imin:imax].squeeze()
	temp = dataset['temperature']['temperature'][kf,:,jmin:jmax,imin:imax].squeeze()
	salt = dataset['salinity']['salinity'][kf,:,jmin:jmax,imin:imax].squeeze()
	u    = dataset['u']['u'][kf,:,jmin:jmax,imin:imax].squeeze()
	v    = dataset['v']['v'][kf,:,jmin:jmax,imin:imax].squeeze()
	spval=zeta.max()
	# write data to netcdf
	currentdate = (ref_roms + dt.timedelta(days=time_vector_corrected[kf])).isoformat()
	ncfile = my_workspace + 'HYCOM_GLBa0.08_' + my_target_grid + '_' + currentdate + '.nc'
	fid = nc.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'HYCOM + NCODA Global 1/12 Analysis (GLBa0.08) \n file created by raphael@esm.rutgers.edu'
        # dimensions
        fid.createDimension('z', nz)
        fid.createDimension('y', ny)
        fid.createDimension('x', nx)
        fid.createDimension('time', None)
        # variables
        latitudes  = fid.createVariable('lat' , 'f8', ('y','x',))
        longitudes = fid.createVariable('lon' , 'f8', ('y','x',))
        zprof      = fid.createVariable('z'   , 'f8', ('z',))
        times      = fid.createVariable('time', 'f8', ('time',))
        o_zeta     = fid.createVariable('zeta', 'f8', ('time','y','x',),fill_value=spval)
        o_temp     = fid.createVariable('temp', 'f8', ('time','z','y','x',),fill_value=spval)
        o_salt     = fid.createVariable('salt', 'f8', ('time','z','y','x',),fill_value=spval)
        o_u        = fid.createVariable('u', 'f8', ('time','z','y','x',),fill_value=spval)
        o_v        = fid.createVariable('v', 'f8', ('time','z','y','x',),fill_value=spval)
        # data
        latitudes[:,:]    = lat_out
        longitudes[:,:]   = lon_out
        zprof[:]          = depth_out
        times[:]          = time_vector_corrected[kf]
        o_zeta[0,:,:]     = zeta
        o_temp[0,:,:,:]   = temp
        o_salt[0,:,:,:]   = salt
        o_u[0,:,:,:]      = u
        o_v[0,:,:,:]      = v
	# attributes
	latitudes.long_name = 'latitude'
	latitudes.units = 'degrees_north'

	longitudes.long_name = 'longitude'
	longitudes.units = 'degrees_east'

	zprof.long_name = 'depth'  
	zprof.units = 'meter'

	times.units = 'days since ' + ref_roms.isoformat().replace('T',' ')
	times.calendar = 'LEAP'

	o_zeta.long_name = 'Sea Surface Height'
	o_zeta.units = 'meters'
	o_zeta.coordinates='lon lat'

	o_temp.long_name = 'Temperature'
	o_temp.units = 'degrees Celcius'
	o_temp.coordinates='lon lat z'

	o_salt.long_name = 'Salinity'
	o_salt.units = 'PSU'
	o_salt.coordinates='lon lat z'

	o_u.long_name = 'Zonal velocity'
	o_u.units = 'm/s'
	o_u.coordinates='lon lat z'

	o_v.long_name = 'Meridional velocity'
	o_v.units = 'm/s'
	o_v.coordinates='lon lat z'
        # close
        fid.close()
	
