#FIX VIP BOUNDARY CONDITIONS

import netCDF4 as nc
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

a = 464
b = 470
# I cord/X COORD
c = 416
d = 425

vip_file  = '/t3/workdir/liz/VIP/VIP_grid/VIP_BRY_y2007.nc' 
v_fid = nc.Dataset(vip_file,'a')
e_temp = v_fid.variables['temp_east']

nums = range(17900,18262+1)
#nums = range(17900,17905)
for nt in nums:

  ct_file = '/data/external/P2/ROMS/CORAL/RUN16/coral_avg_'+str(nt)+'.nc'
  ct_fid = nc.Dataset(ct_file,'r')
  ct_bot_temp = np.squeeze(ct_fid.variables['temp'][:,0,a:b,c:d])
  ct_min = np.min(ct_bot_temp)
  
#  plt.figure()
#  plt.pcolor(np.squeeze(e_temp[nt-nums[0],0:15,100:]))
#  plt.clim(0,25)
#  plt.colorbar()
#  plt.show()
  print(ct_min)

  temp_patch = e_temp[nt-nums[0],0:5,233:271]
  temp_patch = np.where(temp_patch > ct_min,temp_patch,ct_min)
  e_temp[nt-nums[0],0:5,233:271]= temp_patch
#  plt.figure()
#  plt.pcolor(np.squeeze(e_temp[nt-nums[0],0:15,100:])) 
#  plt.clim(0,25)
#  plt.colorbar()
#  plt.show()

v_fid.close()
