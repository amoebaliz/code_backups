import netCDF4 as nc
from matplotlib import cm
import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')
#!ncdump -h Jarvis_5m.grd

#fid = nc.Dataset('/t3/workdir/liz/external_data/JARVIS/BATHYM/jarvis_20m.grd','r')
fid2 = nc.Dataset('/t3/workdir/liz/external_data/JARVIS/BATHYM/jarvis_20m.grd','r')

#x = fid.variables['x'][:]
x2 = fid2.variables['x'][:]
print x2.shape
#y = fid.variables['y'][:]
y2 = fid2.variables['y'][:]

#bathy = fid.variables['z'][:]
bathy2 = fid2.variables['z'][:]
gridx2, gridy2 =np. meshgrid(x2,y2)

#plt.figure()
h = ax.plot_surface(gridx2,gridy2,bathy2, cmap =cm.coolwarm,vmin=-3300,vmax =0, linewidth = 0, shade=False)
fig.colorbar(h, shrink=0.5, aspect=5)
plt.contourf(gridx2,gridy2,bathy2)


#plt.subplot(2,1,1)
#Axes3D.contourf()
#Axes3D.contourf(x,y,bathy,100)
#Axes3D.plot_surface(x,y,bathy,100)
#plt.colorbar()

#plt.subplot(2,1,2)
#ax.contourf(x2,y2,bathy2,100) 
#plt.colorbar()  

plt.show()
