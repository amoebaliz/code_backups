#COMPARING SST FROM VIP and observations
import netCDF4 as nc
import numpy as np
import pyroms
import matplotlib.pylab as plt

grd = pyroms.grid.get_ROMS_grid('VIP')

lonmin_roms = grd.hgrid.lon_vert.min()
lonmax_roms = grd.hgrid.lon_vert.max()
latmin_roms = grd.hgrid.lat_vert.min()
latmax_roms = grd.hgrid.lat_vert.max()

ostia_file = '/t3/workdir/liz/utility_scripts/OSTIA_analysed_sst_2012.nc'
ostia_fid = nc.Dataset(ostia_file)

lon_os = ostia_fid.variables['lon'][:]
lat_os = ostia_fid.variables['lat'][:] 

Ilon = np.where((lon_os > lonmin_roms) & (lon_os < lonmax_roms))
Ilat = np.where((lat_os > latmin_roms) & (lat_os < latmax_roms))

temp_os = ostia_fid.variables['analysed_sst'][:,Ilat[0],Ilon[0]]
temp_os_C = temp_os - 273.15

vip_sst_file = '/t3/workdir/liz/VIP/VIP_runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150205/sst_2007.nc'
vip_fid = nc.Dataset(vip_sst_file)

temp_v = vip_fid.variables['temp'][:]

plt.figure()
plt.subplot(121)
C1 = plt.contourf(np.squeeze(np.mean(temp_os_C,axis = 0)),100)
#plt.clim(27,36)
plt.colorbar(C1)
plt.subplot(122)
C2 = plt.contourf(np.squeeze(np.mean(temp_v,axis = 0)),100)
plt.clim(27,36)
plt.colorbar(C2)
plt.show()

http://data.nodc.noaa.gov/thredds/dodsC/cortad/Version4/cortadv4_row00_col00.nc




