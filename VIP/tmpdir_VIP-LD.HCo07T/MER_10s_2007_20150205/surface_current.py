import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
import pyroms as py

fig = plt.figure(num = 1, figsize = (25,10))
mngr = plt.get_current_fig_manager()
geom = mngr.window.geometry()
x,y,dx,dy = geom.getRect()
mngr.window.setGeometry(200,150,dx,dy)
afreq = 10

ncfile = 'current_avg.nc'
fid = nc.Dataset(ncfile, 'r')
U = fid.variables['u'[:]]
V = fid.variables['v'[:]]
# INTERPOLATE TO RHO-POINTS 
u = np.squeeze((U[:,-1,1:-1,0:-1]+U[:,-1,1:-1,1:])/2)
v = np.squeeze((V[:,-1,0:-1,1:-1]+V[:,-1,1:,1:-1])/2)
u2 = U[:]

u = np.float64(u)
v = np.float64(v)
# CALCULATE CURRENT SPEED 
current_speed = np.sqrt(np.square(u)+np.square(v))
fid.close()

x = np.arange(0,current_speed.shape[1])
y = np.arange(0,current_speed.shape[0])
X, Y = np.meshgrid(x, y)


C = plt.contourf(X,Y,current_speed,100)
plt.clim(0, 0.3)
plt.colorbar()
Q = plt.quiver(X[::afreq,::afreq],Y[::afreq,::afreq],u[::afreq,::afreq],v[::afreq,::afreq],width = 0.0015)



VIP = py.grid.get_ROMS_grid('VIP')
min_row = 70
max_row = 210
rows = np.arange(min_row,max_row+1)
col = 180
X2 = np.tile(rows,(50,1))
U_vals = np.squeeze(u2)
u_slice, u_Z, u_X, u_Y = py.tools.islice(U_vals,col,VIP,'u')


fig = plt.figure()
plt.subplot(131, axisbg = [.5,.5,.5])
cmap = plt.cm.bwr
cmap.set_bad('k',1.)
C1 = plt.contourf(X2,u_Z[:,rows],u_slice[:,rows],100, cmap = cmap)
plt.clim(-.3, 0.3)
plt.ylim(-1000, 0)

#min_row = 125
#max_row = 195
#rows = np.arange(min_row,max_row+1)
col = 280
u_slice, u_Z, u_X, u_Y = py.tools.islice(U_vals,col,VIP,'u')
X2 = np.tile(rows,(50,1))
plt.subplot(132, axisbg = [.5,.5,.5])
C1 = plt.contourf(X2,u_Z[:,rows],u_slice[:,rows],100,cmap=cmap)
plt.clim(-.3, 0.3)
plt.ylim(-1000, 0)

#min_row = 140
#max_row = 210
#rows = np.arange(min_row,max_row+1)
col = 320
u_slice, u_Z, u_X, u_Y = py.tools.islice(U_vals,col,VIP,'u')
X2 = np.tile(rows,(50,1))
plt.subplot(133,axisbg = [.5,.5,.5])
C1 = plt.contourf(X2,u_Z[:,rows],u_slice[:,rows],100,cmap =cmap)
plt.colorbar()
plt.clim(-.3, 0.3)
plt.ylim(-1000, 0) 
plt.show()


