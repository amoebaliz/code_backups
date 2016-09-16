import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import griddata
from datetime import datetime
import pyroms
import pyroms_toolbox
from mpl_toolkits.basemap import Basemap

grd = pyroms.grid.get_ROMS_grid('VIP')

file = '/data/external/P1/Data/Reef_locations/reef_locations_1km.nc'
nc = netCDF4.Dataset(file)
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
lon, lat = np.meshgrid(lon, lat)
coral_cover = nc.variables['coral_cover'][:]

coral_reef = np.zeros(coral_cover.shape)
idx = np.where(coral_cover != 0)
coral_reef[idx] = 1


c =  griddata(lon.flatten(), lat.flatten(), coral_reef.flatten(), grd.hgrid.lon_rho, grd.hgrid.lat_rho)

coral_reef_roms = np.zeros(grd.hgrid.mask_rho.shape)
idx = np.where(c != 0)
coral_reef_roms[idx] = 1
coral_reef_roms = coral_reef_roms * grd.hgrid.mask_rho


# write output in netcdf file
outfile='reef_locations_roms_grd.nc'
nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
nc.Description = 'Coral Reef location mapped on ROMS grid'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'Coral Reef location'

Mp, Lp = grd.hgrid.lon_rho.shape

nc.createDimension('xi_rho', Lp)
nc.createDimension('eta_rho', Mp)

nc.createVariable('lon', 'f8', ('eta_rho', 'xi_rho'))
nc.variables['lon'].long_name = 'longitude'
nc.variables['lon'].units = 'degree_east'
nc.variables['lon'][:] = grd.hgrid.lon_rho

nc.createVariable('lat', 'f8', ('eta_rho', 'xi_rho'))
nc.variables['lat'].long_name = 'latitude'
nc.variables['lat'].units = 'degree_north'
nc.variables['lat'][:] = grd.hgrid.lat_rho

nc.createVariable('coral_reef', 'f4', ('eta_rho', 'xi_rho'), fill_value=1e20)
nc.variables['coral_reef'].long_name = 'coral percentage cover'
nc.variables['coral_reef'].units = 'None'
nc.variables['coral_reef'].valid_min = '0'
nc.variables['coral_reef'].valid_max = '1'
nc.variables['coral_reef'][:] = coral_reef_roms

nc.close()


