import netCDF4 as nc
import numpy as np
import pyroms as py
from matplotlib import pyplot as plt
cmin = 20; cmax = 37; cvals = np.linspace(cmin,cmax,100, endpoint=True)
time = '0279'
row = 177
min_col = 567; max_col = 595
cols = np.arange(min_col,max_col+1)

VIP = py.grid.get_ROMS_grid('VIP')

dir = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo11T/1996/'
data_file = dir + 'VIP-LD.HCo11T_avg_1996-01-01T00:00:00.nc'
fid = nc.Dataset(data_file,'r')
temp = fid.variables['temp'][:]
temp = np.squeeze(temp)

mask_vals = VIP.hgrid.mask_rho
s_rho = VIP.vgrid.s_rho
s_rho = s_rho[:]
mask_grid = mask_vals[:]
X = np.tile(cols,(len(s_rho),1))


# TRANSECT LOCATION PLOT
fig1 = plt.figure()
plt.pcolor(temp[-1,:,:])
plt.plot([min_col,max_col],[row,row],'ko-')
plt.tight_layout(pad=4, w_pad=5, h_pad=1.0)
plt.xlim([0,mask_grid.shape[1]])
plt.ylim([0,mask_grid.shape[0]])
plt.clim(cmin,cmax)
plt.colorbar()
# FIGURE CONFIGURATION
fig2 = plt.figure(num = 2, figsize = (25,7))
ax = fig2.gca()
mngr = plt.get_current_fig_manager()
geom = mngr.window.geometry()
x,y,dx,dy = geom.getRect()
mngr.window.setGeometry(200,50,dx,dy)
# TRANSECT PLOT

temp = fid.variables['temp'][:]
temp = np.squeeze(temp)
temp_slice, temp_Z, temp_X, temp_Y = py.tools.jslice(temp,row,VIP,'rho')

C1 = plt.contourf(X,temp_Z[:,cols],temp_slice[:,cols],200, vmin=cmin,vmax=cmax)
bar_ticks = 11
cbar = fig2.colorbar(C1, ticks=np.linspace(cmin,cmax,bar_ticks,endpoint=True))
cbar.ax.set_ylim([cbar.norm(cmin),cbar.norm(cmax)])
#cbar.outline.set_ydata([cbar.norm(cmin)]*2 + [cbar.norm(cmax)]*4 + [cbar.norm(cmin)]*3)
cbar.ax.set_aspect(40)
plt.show()


# import pyroms
# pyroms.<TAB>
