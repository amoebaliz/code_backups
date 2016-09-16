import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
from mpl_toolkits.basemap import Basemap

slp_dir = '/t3/workdir/liz/external_data/MERRA/MERRA_SLP/'
slp_files = slp_dir + '*'
slp_fid = nc.MFDataset(slp_files,'r')

wind_file = '/t3/workdir/liz/external_data/MERRA/MERRA_WINDS/MERRA_VIP_winds_1979-2013.nc'
wind_fid = nc.Dataset(wind_file,'r')

lats = range(160,240+1)
lons = range(145,234+1)
years = range(1979,2013+1)

ref = dt.datetime(1900,1,1,0,0)

def convert_time(time,ref,incr):
   date = np.zeros(len(time))
   for nt in range(len(time)):
       if incr == 'seconds':
          day_time = ref + dt.timedelta(seconds=time[nt])
       elif incr == 'days':
          day_time = ref + dt.timedelta(days=time[nt])
       date[nt] = pltd.date2num(day_time)
   return date

Uwind = wind_fid.variables['Uwind'][:]
Vwind = wind_fid.variables['Vwind'][:]
wind_time = wind_fid.variables['time'][:]

Uwind = np.squeeze(np.mean(np.mean(Uwind,axis=2),axis=1))
Vwind = np.squeeze(np.mean(np.mean(Vwind,axis=2),axis=1))
wind_plt_time = convert_time(wind_time,ref,'days')

#INDEXING
NEM_index = np.where((Uwind<0) & (Vwind<0))
ENSO_index = np.where((Uwind<0) & (Vwind<0) & (Vwind>-.2))
# MODIFY SLP 
SLP = slp_fid.variables['Pair'][:,lats,lons]
print SLP[NEM_index[0],:,:].shape[0], SLP.shape[0]
composite_vals = np.mean(SLP[NEM_index[0],:,:],axis=0)

print SLP[ENSO_index[0],:,:].shape[0], SLP.shape[0]
composite_vals_2 = np.mean(SLP[ENSO_index[0],:,:],axis=0)
 
plt_lats = slp_fid.variables['lat'][lats]
plt_lons = slp_fid.variables['lon'][lons]

# PLOTTING REGRESSION FIGURE
mlons, mlats = np.meshgrid(plt_lons,plt_lats)
# create new figure, axes instances.
fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])
# setup mercator map projection.
m = Basemap(llcrnrlon=95.,llcrnrlat=-12.,urcrnrlon=160.,urcrnrlat=32.,\
             resolution='l',projection='merc')
im1 = m.pcolormesh(mlons,mlats,composite_vals/100,shading='flat',cmap=plt.cm.jet,latlon=True,vmin=1008,vmax=1014)
m.drawcoastlines()

# draw parallels
m.drawparallels(np.arange(-10,40,10),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
ax.set_title('NEM composite of SLP (hPa), 1979-2007')

cb = m.colorbar(im1,"bottom", size="5%", pad="8%")
#cb.set_label('Regression Slope')

fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])
# setup mercator map projection.
m = Basemap(llcrnrlon=95.,llcrnrlat=-12.,urcrnrlon=160.,urcrnrlat=32.,\
             resolution='l',projection='merc')
im1 = m.pcolormesh(mlons,mlats,composite_vals_2/100,shading='flat',cmap=plt.cm.jet,latlon=True,vmin=1008,vmax=1014)
m.drawcoastlines()

# draw parallels
m.drawparallels(np.arange(-10,40,10),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
ax.set_title('ENSO composite of NEM SLP (hPa), 1979-2007')

cb = m.colorbar(im1,"bottom", size="5%", pad="8%")
#cb.set_label('Regression Slope')


plt.show()
