import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import pyroms as py

# MONTH NAMES
title_val = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

min_col = 190
max_col = 320
cols = np.arange(min_col,max_col)

min_row = 60
max_row = 200
monthly_avg =[]
for J in np.arange(len(title_val)):
     
     ncfile_current = title_val[J] + '_current.nc'
     fid = nc.Dataset(ncfile_current,'r')
     ssh = fid.variables['zeta'[:]]
     ssh = ssh[:]
     ssh_avg = np.squeeze(np.mean(np.mean(ssh[:,min_row:max_row,min_col:max_col],axis=1),axis=0))
     plt.subplot(1,2,1)
     plt.plot(ssh_avg)
     norm_ssh = ssh_avg/(np.amax(ssh_avg))
     monthly_avg =np.append(monthly_avg, norm_ssh, axis=0)
     plt.subplot(1,2,2)
     plt.plot(norm_ssh)
monthly_avg=np.reshape(monthly_avg,(12,len(cols))) 

plt.figure()
C1 = plt.contourf(monthly_avg,100)
plt.show()

# SPATIAL SUBSET

# AVERAGE OVER ALL ROWS

# HOVMOLLER TIME X LON
# matplotlib colored contour
 
