import netCDF4 as nc
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import datetime as dt

yrs = np.arange(1960,2007+1)

months = ['January','February','March','April','May','June',\
'July','August','September','October','November','December']

# FILE DETAILS
file = '/t3/workdir/liz/MODELS/PALAU/palau_aviso_ssh_2010.nc'
fid = nc.Dataset(file,'r')

time_vals = range(60,243)
zeta = np.squeeze(fid.variables['Grid_0001'][time_vals[0],:,:]).T
lats = fid.variables['NbLatitudes'][:]
lons = fid.variables['NbLongitudes'][:]

#zeta = fid.variables['Grid_0001'][:]
#print np.min(zeta)
#print np.max(zeta)
fig = plt.figure(figsize=(9,6))
m = Basemap(llcrnrlon=lons[10],llcrnrlat=lats[20],urcrnrlon=lons[-30],urcrnrlat=lats[-30],resolution='i')
m.drawcoastlines()
m.fillcontinents()
#ax1 = plt.subplot(axisbg=[.2,.2,.2])
#ax1.set_ylim(0,zeta.shape[0])
#ax1.set_xlim(0,zeta.shape[1])
#plt.tight_layout(pad=4.0, w_pad=1, h_pad=4.0)
plt.title('SSH Anomaly \n 1-April-2010 : 31-August-2010')
pcol = m.pcolor(lons, lats[10:-10], zeta[10:-10,:], vmin=-30,vmax=30)
cbar = plt.colorbar(pcol)
cbar.ax.set_ylabel('SLA (cm)', fontsize =14, labelpad = 20, rotation=270)
#ref = dt.datetime(1960,1,1,12,0)
ref_0 = dt.datetime(2010,1,1,0,0)
ref_end = dt.datetime(2010,12,31,0,0)
ims1 = []

# April 1 through end of August 

def get_zeta(i):
#    n = 732+idd
    zeta = np.squeeze(fid.variables['Grid_0001'][i,:,:])
    return zeta.T

#for add in np.arange(5):
for add in time_vals:
#for add in range(50):
    print add
    zeta_val = get_zeta(add)

    ims1.append((m.pcolor(lons,lats[10:-10],zeta_val[10:-10,:],vmin=-30,vmax=30),))

im_ani = animation.ArtistAnimation(fig, ims1, blit=True)
im_ani.save('Palau_2010.gif', writer = 'imagemagick',fps=15)
