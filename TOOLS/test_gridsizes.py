import numpy as np
import netCDF4 as nc

#MaPhil
grid_file = '/t3/workdir/liz/MaPhil/Inputs/Grid/MaPhil_grd_high_res_bathy.nc'

#VIP
#grid_file = '/t3/workdir/liz/VIP/Inputs/Grid/VIP_grd_high_res_bathy_interp.nc'

fid = nc.Dataset(grid_file,'r')
# in the XI direction
pm = fid.variables['pm'][:]
# in the Eta (Y) directon
pn = fid.variables['pn'][:]
# landmask details 
mask_vals = fid.variables['mask_rho'][:]

new_pm = np.ma.masked_where(mask_vals==0,pm)
new_pm = 1/new_pm
print np.ma.max(new_pm)
print np.ma.min(new_pm)
print
new_pn = np.ma.masked_where(mask_vals==0,pn)
new_pn = 1/new_pn
print np.ma.max(new_pn)
print np.ma.min(new_pn)
print




