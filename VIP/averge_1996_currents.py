import netCDF4 as nc
import numpy as np
import pyroms
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    mapid.drawcoastlines(linewidth=0)
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor='0.8', closed=False)
    axs.add_collection(lc)

def get_sst():
    sst = np.mean(np.squeeze(Tfid.variables['temp'][:,:,:,:]),axis=0)
    return sst

def get_vels():
     u_vel = np.mean(np.squeeze(Ufid.variables['u'][:,:,:,:]),axis=0)
     v_vel = np.mean(np.squeeze(Vfid.variables['v'][:,:,:,:]),axis=0)
     # Interpolate
     u_vel_2 = (u_vel[1:-1,1:] + u_vel[1:-1,:-1])/2
     v_vel_2 = (v_vel[1:,1:-1] + v_vel[:-1,1:-1])/2
     # Speed Calculation
     cur_speed = np.sqrt((u_vel_2)**2+(v_vel_2)**2)
     return u_vel_2, v_vel_2, cur_speed

year = 1999
deg_off = 0.02
afreq = 5
cmax = .5
cmin = 0
wx = 121.3
wy = 13.8
VIP = pyroms.grid.get_ROMS_grid('VIP')
#yrs = np.arange(1996,1998+1)

Ilats = [60,232]
Ilons = [168,335]

months = ['JUN','JUL','AUG','SEP']
ndays = [30,31,31,30]

# FILE DETAILS

#Tncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_SST_JJAS_' + str(year) + '.nc'
Tncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_normSST_JJAS_' + str(year) + '.nc'
Uncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_U_JJAS_' + str(year) + '.nc'
Vncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_V_JJAS_' + str(year) + '.nc'  

Tfid = nc.Dataset(Tncfile)
Ufid = nc.Dataset(Uncfile)
Vfid = nc.Dataset(Vncfile)

#vip_file = '/t3/workdir/liz/MODELS/VIP/PostProc/SST_VIP_96-99.nc'
#vip_fid = nc.Dataset(vip_file,'r')
#time_vals = range(60,243)
#sst = np.squeeze(Tfid.variables['temp'][0,:,:,:])
#u = np.squeeze(Ufid.variables['u'][0,:,:,:])
#v = np.squeeze(Vfid.variables['v'][0,:,:,:])
#sst = get_sst()
#u, v,cur_speed = get_vels()
#vip_sst = np.squeeze(vip_fid.variables['temp'][time_vals[0],:,:,:])

#lats = VIP.hgrid.lat_rho[Ilats[0]:Ilats[1],Ilons[0]-5:Ilons[1]+5]
#lons = VIP.hgrid.lon_rho[Ilats[0]:Ilats[1],Ilons[0]-5:Ilons[1]+5]

lats = VIP.hgrid.lat_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
lons = VIP.hgrid.lon_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]

# Figure Details
fig = plt.figure(figsize=(15,7))
ax = plt.subplot(111)
m = Basemap(llcrnrlon=np.min(lons)+.04,llcrnrlat=np.min(lats)+.3,urcrnrlon=np.max(lons)+deg_off,urcrnrlat=np.max(lats)-.19,resolution='h',ax =ax )
#cbar = fig.colorbar(pcol,orientation='horizontal')
u,v,cspeed = get_vels()
pcol = m.pcolor(lons, lats, cspeed, cmap='jet',vmin=cmin,vmax=cmax)
Q = m.quiver(lons[::afreq,::afreq],lats[::afreq,::afreq],u[::afreq,::afreq],v[::afreq,::afreq])
polygon_patch(m,ax)
#tx = plt.text(121.2,13.92,'1999-June-01', fontsize=14)
m.drawmeridians([120.4,121.4])
m.drawparallels([13.4,13.9])
ax.quiver(wx,wy,1,1,headaxislength=3,scale = 2.0,scale_units='inches')

plt.show()
