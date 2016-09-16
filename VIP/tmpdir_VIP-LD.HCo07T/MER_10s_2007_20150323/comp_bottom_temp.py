import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

files = [1,2,250]
title_vals = ['Initial Condition','Day 2', 'Day 250']
lev = np.linspace(0,32,100,endpoint=True)
bgc = [0.7,0.7,0.7]

f, axarr = plt.subplots(3,1, figsize=(12,16), sharex=True)
f.tight_layout(pad=4,w_pad=None,h_pad=None)
f.subplots_adjust(right=0.85,top=0.93)
cbar_ax = f.add_axes([0.9,0.15,0.02,0.7])

for nt in range(len(files)):
  
  file = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150313/VIP-LD.HCo07T_his_' + str(files[nt]).zfill(4) + '.nc'
  fid = nc.Dataset(file,'r')
  bt = fid.variables['temp'][:,0,:,:]
  a = axarr[nt].contourf(np.squeeze(bt),lev)
  axarr[nt].set_axis_bgcolor(bgc)
  axarr[nt].yaxis.tick_left()
  axarr[nt].xaxis.tick_bottom()
  axarr[nt].set_title(title_vals[nt])

plt.suptitle('Bottom Temperatures During 2007 VIP Model Run',fontsize=16)
axarr[1].set_ylabel('I Direction',labelpad=20)
axarr[2].set_xlabel('J Direction',labelpad=20)
plt.colorbar(a,cax=cbar_ax,ticks=[0,8,16,24,32])
cbar_ax.text(-1.9,1.03,'Temperature \n ($^\circ$C)',multialignment='center')
plt.show()
