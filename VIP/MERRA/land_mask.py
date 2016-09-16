import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap




fid = nc.Dataset('/t3/workdir/raphael/from_home/STORAGE/MERRA/mask_MERRA_0-360.nc')

lsm = np.squeeze(fid.variables['lsm'][:])
lat = fid.variables['lat'][:]
lon = fid.variables['lon'][:]

ilon = np.where((lon>117)&(lon<125))
ilat = np.where((lat>10)&(lat<17))

#map = Basemap(projection = 'mill',llcrnrlat=9,urcrnrlat=18,llcrnrlon=116,urcrnrlon=126,resolution='h')


map = Basemap('kav7',lon_0=0,resolution='c')
plt.figure()

map.pcolormesh(lon[ilon[:]],lat[ilat[:]],lsm[np.transpose(ilat[:]),ilon[:]],latlon=True,edgecolor='k')


map.drawcoastlines(linewidth=0.25)

plt.show()

