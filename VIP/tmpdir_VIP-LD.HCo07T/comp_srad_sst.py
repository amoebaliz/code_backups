import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

rad_3hr = '/t3/workdir/liz/VIP/Runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150323/3hr_swrad.nc'
rad_daily = '/t3/workdir/liz/VIP/Runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150401/daily_swrad.nc'

fid_3hr = nc.Dataset(rad_3hr,'r')
fid_daily = nc.Dataset(rad_daily,'r')

sst_3hr = fid_3hr.variables['temp'][:]
sst_daily = fid_daily.variables['temp'][:]

sst_dif = np.squeeze(sst_daily) - np.squeeze(sst_3hr)

f = plt.figure()
ax = f.add_subplot(1,1,1,axisbg = [0.5,0.5,0.5])
p = ax.pcolor(sst_dif,cmap='bwr',vmin=-.1,vmax=0.1)
plt.title('daily - 3hr sw_radiation mean SST (oC)')
plt.xlim(0,sst_dif.shape[1])
plt.ylim(0,sst_dif.shape[0])
plt.colorbar(p)
plt.show()

