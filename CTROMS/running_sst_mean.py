import numpy as np
import netCDF4 as nc
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt

fid = nc.MFDataset('/t3/workdir/liz/MODELS/CTROMS/PostProc/*_CT_SST.nc')

frame_len = 31

ilats = range(400,450)
ilons = range(400,450)

temp = np.squeeze(fid.variables['temp'][:,:,ilats,ilons])
stored_run_max = np.zeros(((temp.shape[0]-frame_len+1),len(ilats),len(ilons)))

for nframe in range(temp.shape[0]-frame_len+1):
    run_max = np.mean(temp[nframe:nframe+frame_len,:,:],axis = 0)
    stored_run_max[nframe,:,:] = np.where(run_max>stored_run_max[nframe,:,:],run_max,stored_run_max[nframe,:,:])

maxm = argrelextrema(np.squeeze(stored_run_max[:,0,0]), np.greater)
print maxm
print (temp.shape[0]-frame_len+1)
mean_val = np.mean(np.squeeze(stored_run_max[:,0,0]))
print mean_val
plt.figure()
plt.plot(np.squeeze(stored_run_max[:,0,0]))
plt.plot([0,(temp.shape[0]-frame_len+1)],[mean_val,mean_val],'r')
plt.show()

