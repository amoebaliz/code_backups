import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from numpy import genfromtxt

my_data=genfromtxt('data.txt', usecols=np.arange(0,2))
my_data = my_data[110:685+1,:]

soi = my_data[:,1]
soi_date = my_data[:,0]
Isoi = np.where(soi<-1)
Isoi2 = np.where(soi>0)
print len(soi)
print len(Isoi[0])
# ZETA DATA
#ilats = np.arange(300,540)
#ilons = np.arange(480,800)

ilats = np.arange(395,465)
ilons = np.arange(630,700)

file = '/t3/workdir/liz/scripts/CT_roms/ctroms_zeta_palau.nc'
fid = nc.Dataset(file,'r')

zeta = fid.variables['zeta'][732:,:,:]
zeta_avg = np.squeeze(np.mean(zeta,axis=0))

zeta_neg_soi = np.mean(zeta[Isoi[0],:,:],axis=0)

#plt.figure()
#plt.plot(zeta_avg)
#plt.plot(soi,'-k')


plt.figure()
plt.pcolor(zeta_avg,vmin=0.3,vmax=0.5)
plt.title('SSH (meters)')
plt.xlim(0,zeta_avg.shape[1])
plt.ylim(0,zeta_avg.shape[0])
plt.colorbar()

#plt.figure()
#plt.pcolor(np.mean(zeta[Isoi2[0],:,:],axis=0),vmin=0.2,vmax=0.9)
#plt.colorbar()
plt.show()
