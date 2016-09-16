import numpy as np 
import netCDF4 as nc
import matplotlib.pyplot as plt

file = '/t1/scratch/liz/Inputs/VIP-LD.HCo11T/Grid/VIP_grd_high_res_bathy_interp2.nc'
fid = nc.Dataset(file)

rho_mask = fid.variables['mask_rho'][:]
angle = fid.variables['angle'][:]

plot_angle = np.ma.masked_where(rho_mask==0,angle)

plt.figure()
plt.pcolor(plot_angle)
plt.xlim(0,plot_angle.shape[1])
plt.ylim(0,plot_angle.shape[0])
plt.colorbar()
plt.show()



