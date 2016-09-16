import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as pltd

west1= [370,492]
west2= [370,465]
west3= [395,465]
west4= [395,492]

east1= [395,451]
east2= [395,480]
east3= [418,480]
east4= [418,451]

#nums = range(18260,18262+1)
nums = range(1,18262+1)
sst_store =np.zeros(len(nums))
ref = dt.datetime(1900,1,1,0,0)
plot_time_store= np.zeros(len(nums))
time_store = np.zeros(len(nums))

for nt in nums:
  file = '/data/external/P2/ROMS/CORAL/RUN16/coral_avg_' + str(nt).zfill(5) + '.nc'
  fid = nc.Dataset(file,'r')
  temp = fid.variables['temp'] 
  time = fid.variables['ocean_time'][:]
  time_store[nt-nums[0]] = time[0]
  newdate = pltd.date2num(ref + dt.timedelta(seconds=time[0]))
  plot_time_store[nt-nums[0]] = newdate
  temp1 = temp[:,-1,west3[1]:west1[1],west1[0]:west3[0]]
  temp1= np.ma.ravel(temp1)
  temp2 = temp[:,-1,east1[1]:east3[1],east1[0]:east3[0]]
  temp2 = np.ma.ravel(temp2)  
 
  vip_sst = np.ma.concatenate((temp2,temp1),axis=1)
  sst_store[nt-nums[0]] = np.mean(vip_sst)

plt.figure()
plt.plot_date(plot_time_store,sst_store,'-k')
plt.show()

fid_new = nc.Dataset('ct_roms_sst_ts.nc','w')
fid_new.createDimension('ocean_time',None)

time = fid_new.createVariable('ocean_time','f8',('ocean_time',))
time.long_name = 'averaged time since initialization'
time.units = 'seconds since 1900-01-01 00:00:00'
time[:] = time_store

plot_time = fid_new.createVariable('plot_time','f8',('ocean_time',))
plot_time[:]=plot_time_store

sst = fid_new.createVariable('ctroms_vip_sst','f8',('ocean_time',))
sst.long_name = 'time-averaged potential temperature'
sst.units = 'Celsius'
sst[:] = sst_store

fid_new.close()


