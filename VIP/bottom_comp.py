import netCDF4 as nc 
import numpy as np
import matplotlib.pyplot as plt

ctroms_file = '/data/external/P2/ROMS/CORAL/RUN16/coral_avg_17901.nc'
ctroms_interp = '/t3/workdir/liz/VIP/VIP_grid/VIP_interpolations/ctroms_temp/coral_avg_17901_VIP.nc'
new_vip_file = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/VIP-LD.HCo07T_his_0001.nc'
old_vip_file = '/t3/workdir/liz/VIP/VIP_runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150304/VIP-LD.HCo07T_his_0001.nc'

levels = np.arange(0,35)

fid1 = nc.Dataset(ctroms_file,'r')
ctroms_temp = fid1.variables['temp'][:,0,450:495,365:435]

fid2 = nc.Dataset(ctroms_interp,'r')
ctinterp_temp = fid2.variables['temp'][:,0,:,:]

fid3 = nc.Dataset(new_vip_file,'r')
new_vip_temp = fid3.variables['temp'][:,0,:,:]

#fid4 = nc.Dataset(old_vip_file,'r')
#old_vip_temp = fid4.variables['temp'][:,0,:,:]

f, axarr = plt.subplots(3,1, figsize=(14,16))
f.tight_layout(pad=4,w_pad=None,h_pad=None)
f.subplots_adjust(right=0.85,top=0.93)
cbar_ax = f.add_axes([0.9,0.15,0.02,0.7])
a = axarr[0].contourf(np.squeeze(ctroms_temp),levels)#; axarr[0,0].set_title('CT-ROMS');axarr[0,0].set_ylabel('Western Boundary')
axarr[1].contourf(np.squeeze(ctinterp_temp),levels)#; axarr[0,1].set_title('VIP')
#axarr[2].contourf(np.squeeze(old_vip_temp),levels, sharex = axarr[1])#; axarr[1,0].set_ylabel('Southern Boundary')
axarr[2].contourf(np.squeeze(new_vip_temp),levels, sharex=axarr[1])

plt.colorbar(a,cax=cbar_ax)
plt.clim(-2,2)
plt.show()

