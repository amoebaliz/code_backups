import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

#bd_fid = nc.Dataset('/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Boundary/hycom_BDRY.nc')
bd_fid = nc.Dataset('/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Boundary/hycom_BDRY.nc')
#bd_fid2 = nc.Dataset('/t3/workdir/liz/MODELS/MAPHIL/Inputs/Boundary/old_boundary/hycom_BDRY.nc')
#temp = np.squeeze(bd_fid.variables['temp_west'][31:31+511,0,:])
temp = np.squeeze(bd_fid.variables['temp_north'][31:31+511,-1,:])
print temp.shape
#temp2 = np.squeeze(bd_fid2.variables['temp_west'][:,0,:])
md_fid = nc.Dataset('/t3/workdir/liz/MODELS/MAPHIL/PostProc/OBCFAC_edit/surf_his_temp.nc')
#md_fid = nc.Dataset('/t3/workdir/liz/MODELS/MAPHIL/Runs/new_inputs/2009/bot_temp.nc')
temp2 = np.squeeze(md_fid.variables['temp'][:,:,-1,:])
print temp2.shape

plt.figure()
plt.pcolor(temp.transpose(),vmin=10,vmax=31)
plt.colorbar()
#plt.ylim(60,130)

plt.figure()
plt.pcolor(temp2.transpose(),vmin=10,vmax=31)
plt.colorbar()
#plt.ylim(60,130)
#plt.xlim(0,31)

plt.show()

