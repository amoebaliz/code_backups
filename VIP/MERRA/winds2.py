# LOOK AT MERRA WINDS OVER VIP
import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as pltd
from mpl_toolkits.basemap import Basemap
import datetime as dt

# ACCESS MERA FILES FOR 2007

u_wind = '/t1/scratch/forcings_sets/MERRA/drowned/drowned_MERRA_Uwind_3hours_2007.nc'
v_wind = '/t1/scratch/forcings_sets/MERRA/drowned/drowned_MERRA_Vwind_3hours_2007.nc'

fig = plt.figure(num = 1, figsize = (25,10))
mngr = plt.get_current_fig_manager()
geom = mngr.window.geometry()
x,y,dx,dy = geom.getRect()
mngr.window.setGeometry(200,150,dx,dy)

# QUIVER ARROW FREQUENCY

Itime = np.arange(0,30*8)
min_lon = 120.5;  max_lon = 122.0
min_lat = 13;   max_lat = 14

fidu = nc.Dataset(u_wind,'r')
fidv = nc.Dataset(v_wind,'r')

lon = fidu.variables['lon'][:]
lat = fidv.variables['lat'][:]
Ilon = np.where((lon > min_lon-.5) & (lon <max_lon+.5))
Ilat = np.where((lat > min_lat-.5) & (lat < max_lat+.5))

Uwind = fidu.variables['Uwind'][:,Ilat[0],Ilon[0]]
Vwind = fidv.variables['Vwind'][:,Ilat[0],Ilon[0]]

fidu.close()
fidv.close()

plot_lon = lon[Ilon[0]]
plot_lat = lat[Ilat[0]]
X, Y = np.meshgrid(plot_lon, plot_lat)

Uwind_avg = np.mean(Uwind[Itime,:,:], axis = 0)
Vwind_avg = np.mean(Vwind[Itime,:,:], axis = 0)

print Uwind_avg.shape
m = Basemap(projection='merc',llcrnrlat=min_lat,urcrnrlat=max_lat,\
            llcrnrlon=min_lon,urcrnrlon=max_lon,lat_ts=20,resolution='h')

m.drawcoastlines()
m.fillcontinents(color='coral',lake_color ='aqua',zorder = 0)
m.drawmapboundary(fill_color='aqua')
Q = m.quiver(X,Y,Uwind_avg,Vwind_avg,latlon = 'true')
        
plt.xlabel('Longitude')
plt.ylabel('Latitude')

plt.show()
