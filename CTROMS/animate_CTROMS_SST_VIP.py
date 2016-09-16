import netCDF4 as nc
import numpy as np
import pyroms
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import datetime as dt

def get_sst(i,fid):
    sst = np.squeeze(fid.variables['temp'][i,:,:,:])
    return sst

CORAL = pyroms.grid.get_ROMS_grid('CORAL')
VIP = pyroms.grid.get_ROMS_grid('VIP')
yrs = np.arange(1996,1998+1)



ct_Ilats = [430,520]
ct_Ilons = [350,430]

time_vals = (np.arange(365)) + 365*2
yrs = np.arange(1996,1998+1)
months = ['January','February','March','April','May','June',\
'July','August','September','October','November','December']


# FILE DETAILS
file = '/t3/workdir/liz/MODELS/CTROMS/PostProc/ct_vip_sst_nointerp.nc'
fid = nc.Dataset(file,'r')
vip_file = '/t3/workdir/liz/MODELS/VIP/PostProc/SST_VIP_96-99.nc'
vip_fid = nc.Dataset(vip_file,'r')
#time_vals = range(60,243)
sst = np.squeeze(fid.variables['temp'][time_vals[0],:,:,:])
vip_sst = np.squeeze(vip_fid.variables['temp'][time_vals[0],:,:,:])

lats = CORAL.hgrid.lat_rho[ct_Ilats[0]:ct_Ilats[1]+1,ct_Ilons[0]:ct_Ilons[1]+1]
lons = CORAL.hgrid.lon_rho[ct_Ilats[0]:ct_Ilats[1]+1,ct_Ilons[0]:ct_Ilons[1]+1]
vip_lats = VIP.hgrid.lat_rho
vip_lons = VIP.hgrid.lon_rho

#print np.min(zeta)
#print np.max(zeta)
fig = plt.figure(figsize=(9,6))
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,0],resolution='h')
m.drawcoastlines()
m.fillcontinents()
#ax1 = plt.subplot(axisbg=[.2,.2,.2])
#ax1.set_ylim(0,zeta.shape[0])
#ax1.set_xlim(0,zeta.shape[1])
#plt.tight_layout(pad=4.0, w_pad=1, h_pad=4.0)
#plt.title('SSH Anomaly \n 1-April-2010 : 31-August-2010')
pcol = m.pcolor(lons, lats, sst, vmin=22,vmax=32)
pcol = m.pcolor(vip_lons, vip_lats, vip_sst, vmin=22,vmax=32)
cbar = plt.colorbar(pcol)
m.plot((vip_lons[0,0],vip_lons[-1,0]), (vip_lats[0,0],vip_lats[-1,0]), '--k')
m.plot((vip_lons[0,0],vip_lons[0,-1]), (vip_lats[0,0],vip_lats[0,-1]), '--k'),\
m.plot((vip_lons[0,-1],vip_lons[-1,-1]), (vip_lats[0,-1],vip_lats[-1,-1]), '--k')
m.plot((vip_lons[-1,-1],vip_lons[-1,0]), (vip_lats[-1,-1],vip_lats[-1,0]), '--k')

#cbar.ax.set_ylabel('SLA (cm)', fontsize =14, labelpad = 20, rotation=270)
#ref = dt.datetime(1960,1,1,12,0)
#ref_0 = dt.datetime(2010,1,1,0,0)
#ref_end = dt.datetime(2010,12,31,0,0)
ims1 = []

for add in time_vals:
#for add in range(5):
    print add
    sst = get_sst(add,fid)
    vip_sst = get_sst(add,vip_fid)
    
    ims1.append((m.pcolor(lons, lats, sst, vmin=22,vmax=32),\
                m.pcolor(vip_lons, vip_lats, vip_sst, vmin=22,vmax=32),))

im_ani = animation.ArtistAnimation(fig, ims1, blit=True)
im_ani.save('CTROMS_sst.gif', writer = 'imagemagick',fps=15)
