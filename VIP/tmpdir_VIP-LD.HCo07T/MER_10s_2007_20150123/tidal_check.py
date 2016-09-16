import netCDF4 as nc
import datetime as dt
import numpy as np
import matplotlib.pylab as plt
import matplotlib

ncfile = 'VIP-LD.HCo07T_sta.nc'
fid = nc.Dataset(ncfile, 'r')
ssh = fid.variables['zeta'[:]]
time = fid.variables['ocean_time'[:]]
temp = fid.variables['temp'[:]]

ssh = ssh[:]
tmp = time[:]
temp = temp[:]
fid.close()

ref = dt.datetime(1900,1,1,0,0)

newdate = []
for t in tmp:
   newdate.append(ref + dt.timedelta(seconds=t))
   





plt.figure()
plt.subplot(211)
plt.plot(newdate,ssh)
plt.subplot(212)
plt.plot(newdate,temp[:,:,-1])

plt.show()
