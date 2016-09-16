import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms

grd = pyroms.grid.get_ROMS_grid('VIP')
depth_val = 50

ncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/1998/VIP-LD.HCo13T_avg_1998-07-24T00:00:00.nc'
fid = nc.Dataset(ncfile)

var = np.squeeze(fid.variables['w'][:,:-1,:,:])

istart = 320
iend = 410
jstart = 190
jend = 260

transect, z, lon, lat = pyroms.tools.transect(var, istart, iend, jstart, jend, grd,Cpos='rho')

print z[:].shape
#w_dept,x,y = pyroms.tools.zslice(var, depth_val, grd, Cpos='w', vert=False, mode='linear')
x = np.tile(np.array(range(transect[:].shape[1])),(transect[:].shape[0],1))

col_val = 0.001

fig = plt.figure()
ax = fig.add_subplot(111,axisbg=[0.5,0.5,0.5])
v = np.linspace(-1*col_val,col_val, 50, endpoint=True)

cs = ax.contourf(x,z,transect,v,extend = 'both', cmap='bwr')
#cs = ax.contourf(x,y,w_dept,v,extend = 'both', cmap='bwr')
cs.cmap.set_over('r')
cs.cmap.set_under('b')
#plt.pcolor(w_dept,vmin=-1*col_val,vmax=col_val,cmap='bwr')
plt.colorbar(cs,ticks=[-.001, 0, 0.001])
plt.show()
