import netCDF4 as nc 
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

fid = nc.Dataset('/t1/scratch/liz/Inputs/VIP-LD.HCo08T/Grid/VIP_grd_high_res_bathy_interp.nc','r')
bathy = fid.variables['h'][:]
mask = fid.variables['mask_rho'][:]

my_data=genfromtxt('/home/liz/ANALYSES/VIP/VIP_stations.txt', usecols=np.arange(0,2),skip_header=2)

x_vals = my_data[:,0]
y_vals = my_data[:,1]

bathy_mask = np.ma.masked_where(mask == 0, bathy)

plt.figure(figsize=(20,10))

plt.pcolor(bathy_mask)
plt.colorbar()

plt.plot(x_vals,y_vals,'o',ms=10, mfc='y',mec='k',mew=2.5)

plt.xlim(0,bathy.shape[1])
plt.ylim(0,bathy.shape[0])

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left='off',        # ticks along the left edge are off
    right='off',       # ticks along the right edge are off
    labelleft='off')   # labels along the left edge are off

plt.show()






