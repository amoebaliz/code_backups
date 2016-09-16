import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


file = '/data/external/P1/Data/Reef_locations/reef_locations_1km.nc'

nc = netCDF4.Dataset(file,'r')
lon = nc.variables['lon_vert'][:] 
lat = nc.variables['lat_vert'][:]
coral_cover = nc.variables['bined_coral_cover'][:]

plt.imshow(coral_cover, interpolation='nearest', origin='lower',extent=[95,175,-25,25])

map = Basemap(projection='merc', llcrnrlon=95,llcrnrlat=-25, urcrnrlon=175, urcrnrlat=25, resolution='f')
lon, lat = np.meshgrid(lon,lat)
x, y = map(lon,lat)

plt.figure(figsize=(12, 9))
map.pcolor(x,y,coral_cover)
plt.show()
