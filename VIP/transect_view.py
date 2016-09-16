import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fid = nc.Dataset('vip_temp_transect.nc','r')
temp = fid.variables['temp'][:]
X = fid.variables['X'][:] 
z = fid.variables['z_field'][:]
plt.figure()

avg_temp = np.squeeze(np.mean(temp,axis=0))
print avg_temp.shape

for nt in range(temp.shape[0]+1):
   plt.pcolor(X,z,np.squeeze(temp[nt,:,:])-avg_temp, cmap='bwr')
   plt.clim(-3,3)
   plt.title(str(nt))
   plt.draw()
   plt.pause(0.1)

#pplot = plt.pcolor(X,z,avg_temp)
#plt.contour(X,z,avg_temp,15,colors='k')
#plt.colorbar(pplot)
#plt.show()
