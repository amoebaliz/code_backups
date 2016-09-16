import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

ncfile = '/t3/workdir/liz/VIP/Runs/VIP-LD.HCo10T/outputs/1996/VIP-LD.HCo10T_his_1996-01-01T12:00:00.nc'

fid = nc.Dataset(ncfile)
u = np.squeeze(fid.variables['u'][:,-1,:,:])
v = np.squeeze(fid.variables['v'][:,-1,:,:])
zeta = np.squeeze(fid.variables['zeta'][:])

v_rho = (v[:-1,:] + v[1:,:])/2
u_rho = (u[:,1:] + u[:,:-1])/2

dv = v_rho[:,2:] - v_rho[:,:-2] 
du = u_rho[2:,:] - u_rho[:-2,:]

fid_grid = nc.Dataset('/t1/scratch/liz/Inputs/VIP-LD.HCo10T/Grid/VIP_grd_high_res_bathy_interp.nc')
x_rho = fid_grid.variables['x_rho'][:]
y_rho = fid_grid.variables['y_rho'][:]
coriolis_f = fid_grid.variables['f'][:]
bathy = fid_grid.variables['h'][:]
dep = bathy+zeta

dx = x_rho[1:-1,2:] - x_rho[1:-1,:-2]  
dy = y_rho[2:,1:-1] - y_rho[:-2,1:-1]

dvdx = dv/dx
dudy = du/dy

dvdx_dudy = dvdx-dudy

pot_vort = (dvdx_dudy)/(coriolis_f[1:-1,1:-1])

fig,ax = plt.subplots()
C = ax.pcolor(pot_vort, cmap=plt.get_cmap('bwr'),vmin=-.1,vmax=.1) 
ax.set_axis_bgcolor((.5, .5, .5))
plt.colorbar(C)
plt.show()

