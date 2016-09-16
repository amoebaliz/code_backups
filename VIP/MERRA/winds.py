import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as pltd
import datetime as dt

# ACCESS MERA FILES FOR 2007

u_wind = '/t1/scratch/forcings_sets/MERRA/drowned/drowned_MERRA_Uwind_3hours_2007.nc'
v_wind = '/t1/scratch/forcings_sets/MERRA/drowned/drowned_MERRA_Vwind_3hours_2007.nc'

fig = plt.figure(figsize = (16,2.5))
# QUIVER ARROW FREQUENCY
afreq = 1

ref = dt.datetime(1900,1,1,0,0)

Itime = np.arange(200*8,300*8)
ufid = nc.Dataset(u_wind,'r')
vfid = nc.Dataset(v_wind,'r')

min_lon = 120.5;  max_lon = 122.0
min_lat = 13;   max_lat = 14

hr_time = ufid.variables['time'][:]
lon = ufid.variables['lon'][:]
lat = vfid.variables['lat'][:]

Ilon = np.where((lon > min_lon) & (lon <max_lon))
Ilat = np.where((lat > min_lat) & (lat < max_lat))

Uwind = ufid.variables['Uwind'][:,Ilat[0],Ilon[0]]
Vwind = vfid.variables['Vwind'][:,Ilat[0],Ilon[0]]

day_u = np.zeros((365,Uwind.shape[1],Uwind.shape[2]))
day_v = np.zeros((365,Uwind.shape[1],Uwind.shape[2]))
date = np.zeros(365)
for nt in range(0,365):
   day_u[nt,:,:] = np.mean(Uwind[8*(nt):8*(nt+1),:,:],axis=0)
   day_v[nt,:,:] = np.mean(Vwind[8*(nt):8*(nt+1),:,:],axis=0)
   date_val = np.asscalar(np.mean(hr_time[8*nt:8*(nt+1)]))
   day_time = ref + dt.timedelta(days=date_val)
   date[nt] = pltd.date2num(day_time) 

avg_u = np.mean(np.mean(day_u,axis=2),axis=1)
avg_v = np.mean(np.mean(day_v,axis=2),axis=1)

Y2 = np.zeros(len(avg_u))
#plt.figure()
Q = plt.quiver(date[::afreq],Y2[::afreq],avg_u[::afreq],avg_v[::afreq],width=0.001, scale=90)
qk = plt.quiverkey(Q,.05,.84,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})
plt.gca().xaxis_date()
plt.gca().axes.get_xaxis().tick_bottom()
plt.gca().axes.get_yaxis().set_visible(False)
box_hh = 0.6
plt.plot([date[109],date[109]],[-1*box_hh,box_hh],'-r')
plt.plot([date[183],date[183]],[-1*box_hh,box_hh],'-r')
plt.plot([date[109],date[183]],[box_hh,box_hh],'-r')
plt.plot([date[109],date[183]],[-1*box_hh,-1*box_hh],'-r')

plt.plot([date[214],date[214]],[-1*box_hh,box_hh],'-b')
plt.plot([date[230],date[230]],[-1*box_hh,box_hh],'-b')
plt.plot([date[214],date[230]],[box_hh,box_hh],'-b')
plt.plot([date[214],date[230]],[-1*box_hh,-1*box_hh],'-b')

plt.plot([date[274],date[274]],[-1*box_hh,box_hh],'-b')
plt.plot([date[289],date[289]],[-1*box_hh,box_hh],'-b')
plt.plot([date[274],date[289]],[box_hh,box_hh],'-b')
plt.plot([date[274],date[289]],[-1*box_hh,-1*box_hh],'-b')

plt.xlim(date[0]-30,date[-1]+10)
plt.ylim(-.7,.7)
plt.subplots_adjust(bottom=0.15)
plt.tight_layout(pad = 2.0)

plt.title("Daily Average 2m Wind Speed/Direction Over the VIP (MERRA, m/s)")


















plt.show()
