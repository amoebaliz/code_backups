import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from mpl_toolkits.basemap import Basemap
import pyroms

#file = 'CORAL_hires/coral_cover_CORALx4.nc'
file = '/data/external/P1/Data/Reef_locations/reef_locations_roms_grd.nc'

nc = netCDF4.Dataset(file,'r')
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
coral_cover = nc.variables['coral_reef'][:]
nc.close()

idx = np.where(coral_cover != 0)
reef = np.zeros(coral_cover.shape)
reef[idx] = 1
reef = np.ma.masked_where(reef == 0, reef)

grd = pyroms.grid.get_ROMS_grid('VIP')

map = pyroms.utility.get_grid_proj(grd)
x, y = map(lon,lat)

plt.figure(figsize=(12, 9))

map.drawcoastlines()
map.drawmeridians([100,120,140,160], labels=[0,0,0,1], fmt='%d', fontsize=18)
map.drawparallels([-20,-10,0,10,20], labels=[1,0,0,0], fmt='%d', fontsize=18)
map.etopo(alpha=0.8)

#Basemap.plot(map,x,y,coral_cover,facecolors='Magenta',edgecolors='Magenta')
#Basemap.pcolor(map,x,y,coral_cover,facecolors='Magenta',edgecolors='Magenta')
Basemap.scatter(map,x[idx],y[idx], s=5, marker='o', facecolors='Magenta',edgecolors='Magenta')

xv, yv = map(grd.hgrid.lon_vert, grd.hgrid.lat_vert)
plt.plot(xv[0,:],yv[0,:],ls='-',lw=1.5,color='OrangeRed')
plt.plot(xv[-1,:],yv[-1,:],ls='-',lw=1.5,color='OrangeRed')
plt.plot(xv[:,0],yv[:,0],ls='-',lw=1.5,color='OrangeRed')
plt.plot(xv[:,-1],yv[:,-1],ls='-',lw=1.5,color='OrangeRed')

plt.show()

#outfile='reef_locations.png'
#plt.savefig(outfile, dpi=150, facecolor='w', edgecolor='w', \
#                        orientation='portrait')

