# 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt 
from scipy import stats



nt = 5

fid = nc.Dataset('ct_roms_sst_ts.nc','r')
sst = fid.variables['ctroms_vip_sst'][:]
print sst.shape
#sst_store = np.zeros(len(sst)-nt+1)

slope, intercept, r_value, p_value, std_err = stats.linregress(range(0,len(sst)),sst)

print format(p_value,'.10f')
print format(10*365*slope, 'f')

#for n in range(nt/2,len(sst)-2):
# nums = range(nt/2,10)
# for n in nums:
#  print n 
#  sst_store[n-nums[0]] = np.mean(sst[n-2:n+2])
 



