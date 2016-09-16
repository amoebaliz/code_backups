import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as pltd
import datetime as dt
import pyroms as py
import pyroms.tools as pyt

ref = dt.datetime(1900,1,1,0,0)
#min_lon = 120;  max_lon = 122.0
#min_lat = 13;   max_lat = 14.5
#Expanded Domain
min_lon = 117;  max_lon = 125
min_lat = 10;   max_lat = 17

Uwind = []
Vwind = []
date = []
for nt in range(1996,1999+1):
    print nt
    # ACCESS MERRA FILES FOR 1996-1999
    u_wind = '/t1/scratch/forcings_sets/MERRA/drowned/drowned_MERRA_Uwind_3hours_' + str(nt) + '.nc'
    v_wind = '/t1/scratch/forcings_sets/MERRA/drowned/drowned_MERRA_Vwind_3hours_' + str(nt) + '.nc'
    ufid = nc.Dataset(u_wind,'r')
    vfid = nc.Dataset(v_wind,'r')

    sub_hr_time = ufid.variables['time'][:]
    lon = ufid.variables['lon'][:]
    lat = vfid.variables['lat'][:]

    Ilon = np.where((lon > min_lon) & (lon <max_lon))
    Ilat = np.where((lat > min_lat) & (lat < max_lat)) 
    uwind = ufid.variables['Uwind'][:,Ilat[0],Ilon[0]]
    vwind = vfid.variables['Vwind'][:,Ilat[0],Ilon[0]]
   
    nday = len(sub_hr_time)/8
    day_u = np.zeros((nday,uwind.shape[1],uwind.shape[2]))
    day_v = np.zeros((nday,uwind.shape[1],uwind.shape[2]))
    day_date = np.zeros(nday)

    for jt in range(nday):
        day_u[jt,:,:] = np.mean(uwind[8*(jt):8*(jt+1),:,:],axis=0)
        day_v[jt,:,:] = np.mean(vwind[8*(jt):8*(jt+1),:,:],axis=0)
        date_val = np.asscalar(np.nanmean(sub_hr_time[8*jt:8*(jt+1)]))
        day_time = ref + dt.timedelta(days=date_val) 
        day_date[jt] = pltd.date2num(day_time)

    avg_u = np.mean(np.mean(day_u,axis=2),axis=1)
    avg_v = np.mean(np.mean(day_v,axis=2),axis=1)        

    Uwind = np.append(Uwind,avg_u)
    Vwind = np.append(Vwind,avg_v)
    date = np.append(date, day_date)
# QUIVER ARROW FREQUENCY
afreq = 1
Y2 = np.zeros(len(date))
# calculate surface transport

# TRANSECT POINTS [i,j]
# vip1 = [278,144]; vip2 = [278,178]
# DEPTH RANGES
# h1 = 0; h2 = -50

# VIP = py.grid.get_ROMS_grid('VIP')
# vip_dir = '/t3/workdir/liz/VIP/Runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150401/'

#vip_nums = range(5,361+1)
vip_dir = '/t3/workdir/liz/scripts/VIP_analyses/vip_u_transport_col278_rows144-178.nc'
vip_fid = nc.Dataset(vip_dir)


#def get_transports(i):
#    ncfile = vip_dir + 'VIP-LD.HCo07T_avg_' +str(vip_nums[i]).zfill(4) + '.nc'
#    fid = nc.Dataset(ncfile,'r')
#    u = fid.variables['u'][:]; u.setflags(write=False)
#    v = fid.variables['v'][:]        
#    time = fid.variables['ocean_time'][:]
    
#    v = np.squeeze(v)
#    posu = np.squeeze(np.ma.copy(u)) 
#    negu = np.squeeze(np.ma.copy(u))   
#    posu[posu<0]=0 
#    negu[negu>0]=0
#    day_time = ref + dt.timedelta(seconds=time[0])
    
#    pos_ut,transpv = pyt.section_transport_z(posu, v, VIP, vip1[0], vip2[0], vip1[1], vip2[1], h1, h2)
#    neg_ut,transpv = pyt.section_transport_z(negu, v, VIP, vip1[0], vip2[0], vip1[1], vip2[1], h1, h2)
#    date2 = pltd.date2num(day_time)
#    return pos_ut, neg_ut, time, date2

# Building Plot values 
#pos_u_trans = np.zeros(len(vip_nums))
#neg_u_trans = np.zeros(len(vip_nums))
#date2 = np.zeros(len(vip_nums))
#delta_time = np.zeros(len(vip_nums))

pos_u_trans = vip_fid.variables['pos_transport'][:-1]
neg_u_trans = vip_fid.variables['neg_transport'][:-1]
time  = vip_fid.variables['ocean_time'][:-1]
date2 = np.zeros(len(time))
for nt in range(len(time)):
   vip_time = ref + dt.timedelta(days=time[nt])
   date2[nt] = pltd.date2num(vip_time)
#for nt in range(len(vip_nums)):
#    print nt
#    pos_u_trans[nt], neg_u_trans[nt], delta_time[nt], date2[nt]  = get_transports(nt)[:]

#print 'save to file'

# SAVE TO NETCDF file
#file = 'vip_u_transport_col278_rows144-178.nc'
#ncid = nc.Dataset(file,'w')
#ncid.createDimension('ocean_time', None)
#ncid.title = 'zonal transports between ' + str(h1) + ':' + str(-1*h2) + ' meters, at VIP column 278 and rows 144:178'

#ocean_time = ncid.createVariable('ocean_time', 'f8',('ocean_time',))
#ocean_time.units = 'seconds since 1900-01-01 00:00:00'
#ocean_time[:] = delta_time[:]

#pos_transport = ncid.createVariable('pos_transport','f8',('ocean_time',))
#pos_transport.units = 'm3/s'
#pos_transport.long_name = 'eastward transport'
#pos_transport[:] = pos_u_trans[:]

#neg_transport = ncid.createVariable('neg_transport','f8',('ocean_time',))
#neg_transport.units = 'm3/s'
#neg_transport.long_name = 'westward transport'
#neg_transport[:] = neg_u_trans[:]
#ncid.close()

### FIGURE GENERATION
fig,ax = plt.subplots(1, figsize = (16,5))
fSv = 1000000
# QUIVER ARROW FREQUENCY
#ax[0].axes.get_xaxis().tick_bottom()
ax.plot(date2,pos_u_trans/fSv,'-b')
ax.plot(date2,neg_u_trans/fSv,'-r')
ax.set_xlim(date2[0],date2[-1])
ax.xaxis_date()
ax.set_ylim(-.45,.45)
#plt.subplots_adjust(bottom=0.15)
#plt.tight_layout(pad = 2.0)
#plt.title("Daily Average 2m Wind Speed/Direction Over the VIP (MERRA, m/s)")
afreq = 1
ax2 = ax.twinx()
Q = ax2.quiver(date[::afreq],Y2[::afreq],Uwind[::afreq],Vwind[::afreq],width=0.001, scale=100)
qk = ax2.quiverkey(Q,.05,.84,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})
ax2.xaxis_date()
ax2.axes.get_yaxis().set_visible(False)
plt.title("Eastward (blue) and Westward (red) Transports in the upper 50m (Sv)")

plt.show()
