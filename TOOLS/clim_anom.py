import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


# Calculate climatology and 
#

fid = nc.Dataset('ctroms_u_transport_col395_monthly.nc','r')

pos_T = fid.variables['pos_transport'][:]
neg_T = fid.variables['neg_transport'][:]

pos_T2 = pos_T.reshape(-1,12)
neg_T2 = neg_T.reshape(-1,12)

tot_T2 = pos_T2+neg_T2

clim = np.zeros(12)
stdev  = np.zeros(12)
for nt in range(12):
    clim[nt] = np.mean(tot_T2[:,nt],axis=0)
    stdev[nt] = np.std(tot_T2[:,nt],axis=0) 

fSv = 1000000
plt.figure()
plt.plot(clim/fSv)
plt.plot((clim-stdev)/fSv,'--k')
plt.plot((clim+stdev)/fSv,'--k')
plt.show()
# if not 1-D array, unravel it. 

