import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

fid = nc.Dataset('/t1/scratch/liz/tmpdir_VIP-LD.HCo13T/outputs/1996/VIP_sst.nc')
temp = fid.variables['temp'][:]


for nt in range(temp.shape[0]):
    avg_temp = np.mean(temp[nt,:,:,:])
    plt.figure()
    plt.pcolor(np.squeeze(temp[nt,:,:,:])-avg_temp,vmin=-1,vmax=1)
    plt.colorbar()
    plt.show()

