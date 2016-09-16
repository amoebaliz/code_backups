# calculating transports in the VIP 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms.tools as pyt
import pyroms as py
import datetime as dt
import matplotlib.dates as pltd

vip_dir = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo11T'
# vip_dir = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/'
#ctroms_dir = '/data/external/P2/ROMS/CORAL/RUN16/'
ctroms_interp_dir = ''
ref = dt.datetime(1900,1,1,0,0)
# TRANSECT POINTS [i,j]
vip1 = [170,89]; vip2 = [186,66]

# GRIDS
VIP = py.grid.get_ROMS_grid('VIP')
CORAL = py.grid.get_ROMS_grid('CORAL')

#vip_nums = np.arange(5,20 + 1)
vip_nums = np.arange(5, 217 + 1)
transect_store = np.ma.zeros((len(vip_nums),50,24))
vip_transport = np.zeros(len(vip_nums))
time1 = np.zeros(len(vip_nums))
time = np.zeros(len(vip_nums))
n = 0
for nt in vip_nums:

   vip_file = vip_dir + 'VIP-LD.HCo07T_avg_' + str(nt).zfill(4) + '.nc'
   vip_fid = nc.Dataset(vip_file,'r')
   # u = vip_fid.variables['u'][:]; u = np.squeeze(u)   
   # v = vip_fid.variables['v'][:]; v = np.squeeze(v)
   temp = vip_fid.variables['temp'][:]; temp = np.squeeze(temp)
   time[nt-vip_nums[0]] = vip_fid.variables['ocean_time'][:]
   s_rho = vip_fid.variables['s_rho'][:];vs_rho = s_rho[:]

   # u_up = np.where( 
   # u_low =
   # transpu, transpv = pyt.section_transport(u, v, VIP,vip1[0],vip2[0],vip1[1],vip2[1]) 
   # time1[nt-vip_nums[0]] = pltd.date2num(ref + dt.timedelta(seconds=time[0])) 
   # vip_transport[nt-vip_nums[0]] = transpu + transpv
   # print transpu+transpv
   transect_store[nt-vip_nums[0],:,:], z, lon, lat = pyt.transect(temp,vip1[0],vip2[0],vip1[1],vip2[1],VIP,'rho',spval=1e+37)
   
fid2 = nc.Dataset('vip_temp_transect.nc','w')
#print z.shape
fid2.createDimension('ocean_time', None)
fid2.createDimension('depth', 50)
fid2.createDimension('lateral', 24)


#p1 = fid2.createVariable('p1','f8')
#p1[:]  = vip1
#p2 = fid2.createVariable('p2','f8')
#p2[:]  = vip2
X = np.tile(np.arange(0,z.shape[1]),(len(s_rho),1))
X_var = fid2.createVariable('X','f8',('depth','lateral',))
X_var[:] = X
ocean_time = fid2.createVariable('ocean_time','f8',('ocean_time',))
ocean_time.long_name = 'averaged time since initialization'
ocean_time.units = 'seconds since 1900-01-01 00:00:00'
ocean_time[:] = time
z_field = fid2.createVariable('z_field','f8',('depth','lateral',))
z_field.units = 'meter'
z_field[:] = z
temp = fid2.createVariable('temp','f8',('ocean_time','depth','lateral',))
temp.long_name = 'time-averaged potential temperature' 
temp.units = 'Celsius'
print transect_store.shape
temp[:] = transect_store
fid2.close()

X = np.tile(np.arange(0,z.shape[1]),(len(s_rho),1))
plt.figure()
plt.pcolor(X,z,np.squeeze(np.mean(transect_store,axis=0)))
plt.colorbar()
plt.show()


ct_nums = np.arange(17900,18262+1)
#ct_nums = np.arange(17900,17901+1)
ctroms_transport = np.zeros(len(ct_nums))
transect_store = np.ma.zeros((len(ct_nums),50,5))

time2 = np.zeros(len(ct_nums))
for nt in ct_nums:
   ctroms_file = ctroms_dir + 'coral_avg_' + str(nt) + '.nc'
   ct_fid = nc.Dataset(ctroms_file,'r')
   u = ct_fid.variables['u'][:]; u = np.squeeze(u)
   v = ct_fid.variables['v'][:]; v = np.squeeze(v)
   s_rho = ct_fid.variables['s_rho'][:]; v = np.squeeze(v)
   time = ct_fid.variables['ocean_time'][:]; time=time[:]
   
   transpu, transpv = pyt.section_transport(u, v, CORAL,ct1[0],ct2[0],ct1[1],ct2[1])
   ctroms_transport[nt-ct_nums[0]] = transpu + transpv
   time2[nt-ct_nums[0]] = pltd.date2num(ref + dt.timedelta(seconds=time[0]))
fid2.createVariable('z','f8',('depth','lateral',))
fid2.createVariable('z','f8',('depth','lateral',))
fid2.createVariable('z','f8',('depth','lateral',))
fid2.createVariable('z','f8',('depth','lateral',))
   #print transpu+transpv
   #transect_store[n,:,:], z, lon, lat = pyt.transect(u,ct1[0],ct2[0],ct1[1],ct2[1],CORAL,'u',spval=1e+37)
   #n = n+1

X = np.tile(np.arange(0,z.shape[1]),(len(s_rho),1))
plt.figure()
plt.pcolor(X,z,np.squeeze(np.mean(transect_store,axis=0)))
plt.colorbar()

fig = plt.figure() 
plt.plot_date(time2,ctroms_transport,'-k')
plt.plot_date(time1,vip_transport,'-b')
fig.autofmt_xdate()
plt.show()
