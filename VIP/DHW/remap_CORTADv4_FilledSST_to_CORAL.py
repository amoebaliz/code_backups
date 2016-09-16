import numpy as np
import netCDF4 as netCDF
from datetime import datetime

import pyroms
import pyroms_toolbox


print 'Load CORTADv4 file'
nc_data = netCDF.Dataset('/data/external/P1/Data/CORTAD/Version4/cortadv4_row03_col13.nc')
data = nc_data.variables['sst']
time = nc_data.variables['time'][:]

start_time = pyroms_toolbox.date2jday(datetime(1981,01,01))
end_time = pyroms_toolbox.date2jday(datetime(2012,01,01))
tidx = np.where((time >= start_time) & (time <= end_time))

# load coral grid object
grd = pyroms.grid.get_ROMS_grid('VIP')

# define some variables
wts_file = 'remap_weights_CORTADv4_FilledSST_to_CORAL_bilinear.nc'
nt = len(tidx[0])
Mp, Lp = grd.hgrid.mask_rho.shape
sstr = np.zeros((Mp,Lp))

# create runoff file
file = 'remapped_VIP_CORTADv4_FilledSST.nc'
nc = netCDF.Dataset(file, 'w', format='NETCDF3_64BIT')
nc.Description = 'CORTADv4 filled SST data remapped on ROMS VIP grid'
nc.Author = 'remap_CORTADv4_to_VIP.py'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'CORTADv4 SST data'

# creat dimensions and variables
nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
nc.createDimension('time', None)

nc.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
nc.variables['lon_rho'].long_name = 'longitude of RHO-points'
nc.variables['lon_rho'].units = 'degree_east'
nc.variables['lon_rho'].field = 'lon_rho, scalar'
nc.variables['lon_rho'][:] = grd.hgrid.lon_rho

nc.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
nc.variables['lat_rho'].long_name = 'latitude of RHO-points'
nc.variables['lat_rho'].units = 'degree_north'
nc.variables['lat_rho'].field = 'lat_rho, scalar'
nc.variables['lat_rho'][:] = grd.hgrid.lat_rho

nc.createVariable('time', 'f8', ('time'))
nc.variables['time'].long_name = 'time'
nc.variables['time'].units = 'days since 1900-01-01 00:00:00'
nc.variables['time'][:] = time[tidx]

nc.createVariable('sst', 'f8', ('time', 'eta_rho', 'xi_rho'))
nc.variables['sst'].units = 'Celsius'
nc.variables['sst'].long_name = 'CORTADv4 SST'
nc.variables['sst'].land_value = -1000
nc.variables['sst'].missing_pixel_value = -10000


mask_idx = np.where(grd.hgrid.mask_rho == 0)

nct=0
for t in tidx[0]:
    print 'Remapping SST for time %f' %time[nct]
    # bilinear horizontal interpolation using scrip
    sstc = data[t,:,:]
    #cloud_idx = np.where(sstc == -10000)
    #sstc[cloud_idx] = -1e30
    sstr = pyroms.remapping.remap(sstc, wts_file, \
                                           spval=-1000)

    # cloud
    #cloud_idx = np.where(sstr <= 0)
    #sstr[cloud_idx] = -10000

    # spval
    sstr[mask_idx] = -1000

    # write data in destination file
    nc.variables['sst'][nct] = sstr

    nct = nct + 1



# close netcdf file
nc.close()
                                                     
