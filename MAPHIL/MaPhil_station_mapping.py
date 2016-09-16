import netCDF4 as nc 
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

fid = nc.Dataset('/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Grid/MaPhil_grd_high_res_bathy_mixedJerlov.nc')
bathy = fid.variables['h'][:]
mask = fid.variables['mask_rho'][:]
lats = fid.variables['lat_rho'][:]
lons = fid.variables['lon_rho'][:]

#my_data=genfromtxt('MaPhil_stations.txt', usecols=np.arange(0,2),skip_header=2)

my_data =   np.array([[124.6293029,   10.91280546], \
             [124.6608268,   10.88025355],\
             [124.6820712,   10.83845005],\
             [124.7142805,   10.80315692],\
             [124.7372381,   10.76409463],\
             [124.7478603,   10.71920727],\
             [124.7468324,   10.67294929],\
             [124.7266159,   10.63183109],\
             [124.7197629,   10.59002759],\
             [124.7125672,   10.54376961],\
             [124.6200513,   10.8076114],\
             [124.6358132,   10.67603316],\
             [124.5693388,   10.56364341],\
             [124.786665,    10.743685]])
       
x_vals = my_data[:,0]
y_vals = my_data[:,1]

bathy_mask = np.ma.masked_where(mask == 0, bathy)

plt.figure()

plt.pcolor(lons,lats,bathy_mask)
plt.colorbar()

plt.plot(x_vals,y_vals,'o',mfc='y',mec='k',mew=1.5)
plt.plot([124.785,124.784,124.786],[10.7431,10.7436,10.7427],'o', mfc='g',mec='k',mew=1.5)
plt.xlim(np.min(lons),np.max(lons))
plt.ylim(np.min(lats),np.max(lats))

plt.show()






