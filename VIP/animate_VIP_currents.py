import netCDF4 as nc
import numpy as np
import pyroms
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
    #mapid.drawmapboundary(fill_color=[.9,.97,1])
    mapid.drawcoastlines(linewidth=0)
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor='0.8', closed=False)
    axs.add_collection(lc)

#def get_sst():
def get_sst(i):
    sst = np.squeeze(Tfid.variables['temp'][i,:,:,:])
    #sst = np.mean(np.squeeze(Tfid.variables['temp'][:]),axis=0)
    return sst

def get_vels2():
     u_vel = np.mean(np.squeeze(Ufid.variables['u'][:]),axis=0)
     v_vel = np.mean(np.squeeze(Vfid.variables['v'][:]),axis=0)

     u_vel_2 = (u_vel[1:-1,:-1]+u_vel[1:-1,1:])/2
     v_vel_2 = (v_vel[:-1,1:-1]+v_vel[1:,1:-1])/2

     u4 = u_vel_2*(np.cos(angs[1:-1,1:-1])) + v_vel_2*(np.cos(np.pi/2 + angs[1:-1,1:-1]))
     v4 = u_vel_2*(np.sin(angs[1:-1,1:-1])) + v_vel_2*(np.sin(np.pi/2 + angs[1:-1,1:-1]))
     print u4.shape, v4.shape
     # Speed Calculation
     cur_speed = np.sqrt((u_vel_2)**2+(v_vel_2)**2)

     return u4, v4, cur_speed
def get_vels(i):
     u_vel = np.squeeze(Ufid.variables['u'][i,:,:,:])
     v_vel = np.squeeze(Vfid.variables['v'][i,:,:,:])
     # Interpolate
     u_vel_2 = (u_vel[1:-1,:-1]+u_vel[1:-1,1:])/2
     v_vel_2 = (v_vel[:-1,1:-1]+v_vel[1:,1:-1])/2
     # projection
     u4 = u_vel_2*(np.cos(angs[1:-1,1:-1])) + v_vel_2*(np.cos(np.pi/2 + angs[1:-1,1:-1]))
     v4 = u_vel_2*(np.sin(angs[1:-1,1:-1])) + v_vel_2*(np.sin(np.pi/2 + angs[1:-1,1:-1]))

     # Speed Calculation
     cur_speed = np.sqrt((u_vel_2)**2+(v_vel_2)**2)
     print u4.shape, v4.shape, cur_speed.shape

     return u4, v4, cur_speed

def get_wind(i):
    u_w = np.mean(WUfid.variables['Uwind'][i,:,:])
    v_w = np.mean(WVfid.variables['Vwind'][i,:,:])
    return u_w, v_w
MEEP = 2
year = 1999
deg_off = 0.02
afreq = 5
cmax = .8
#cmax = 
cmin = 0
wx = 121.3
wy = 13.8
VIP = pyroms.grid.get_ROMS_grid('VIP')

Ilats = [60,232]
Ilons = [168,335]

months = ['JUN','JUL','AUG','SEP']
ndays = [30,31,31,30]

# FILE DETAILS

#Tncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_SST_JJAS_' + str(year) + '.nc'
#Tncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_normSST_JJAS_' + str(year) + '.nc'
Uncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_U_JJAS_' + str(year) + '.nc'
Vncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(year) + '/VIP_V_JJAS_' + str(year) + '.nc'
WUfile = '/t3/workdir/liz/external_data/MERRA/Uwind_MERRA_daily_' + str(year) + '.nc'
WVfile = '/t3/workdir/liz/external_data/MERRA/Vwind_MERRA_daily_' + str(year) + '.nc'

#Tfid = nc.Dataset(Tncfile)
Ufid = nc.Dataset(Uncfile)
Vfid = nc.Dataset(Vncfile)
WUfid = nc.Dataset(WUfile)
WVfid = nc.Dataset(WVfile)
 
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

lats = VIP.hgrid.lat_rho#[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
lons = VIP.hgrid.lon_rho#[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
angs = VIP.hgrid.angle_rho

print lats.shape,lons.shape

lats2 = lats[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]                     
lons2 = lons[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
# Figure Details
fig = plt.figure(figsize=(15,7))
ax = plt.subplot(111)
m = Basemap(llcrnrlon=np.min(lons2)+.04,llcrnrlat=np.min(lats2)+.3,urcrnrlon=np.max(lons2)+deg_off,urcrnrlat=np.max(lats2)-.19,resolution='h',ax =ax )
#cbar = fig.colorbar(pcol,orientation='horizontal')
if MEEP == 1:
   #sst = get_sst(0)
   u,v,cur_speed = get_vels(0)

   #sst2  = sst[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
   #sst2[sst2>0]=0
   # shift one down and to left because of "trimming"
   u2    = u[Ilats[0]-1:Ilats[1],Ilons[0]-1:Ilons[1]]
   v2    = v[Ilats[0]-1:Ilats[1],Ilons[0]-1:Ilons[1]]
   cs2   = cur_speed[Ilats[0]-1:Ilats[1],Ilons[0]-1:Ilons[1]] 

   pcol = m.pcolor(lons2, lats2, cs2, cmap='jet',vmin=cmin,vmax=cmax)
   #pcol = m.pcolor(lons2, lats2, sst2, cmap='bwr',vmin=cmin,vmax=cmax)
   Q = m.quiver(lons2[::afreq,::afreq],lats2[::afreq,::afreq],u2[::afreq,::afreq],v2[::afreq,::afreq])
   polygon_patch(m,ax)
   tx = plt.text(121.2,13.92,'1999-June-01', fontsize=14)
   m.drawmeridians([120.4,121.4])
   m.drawparallels([13.4,13.9])

imon = 0
mon_day = 1
# ANIMATION
def updatefig(i):
    global im1,im2,im3,im4,tx,imon,mon_day

    if i>ndays[0]-1:
       mon_day = i+1 - np.sum(ndays[:imon])
    else:
       mon_day = i+1
    if mon_day == ndays[imon]+1:
       if imon == len(months)-1:
          imon = 0
       else:
          imon += 1
       mon_day = 1
    mon = months[imon]    
    day = mon_day      

    # REMOVE images after first step
    if i > 0:
       im1.remove()
       im2.remove()
       im4.remove()
    #sst = get_sst(i)
    u, v,cur_speed = get_vels(i)
    u_w, v_w = get_wind(i+152)

    #sst2  = sst[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
    #sst2[sst2>0]=0
    u2    = u[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
    v2    = v[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1] 
    cs2   = cur_speed[Ilats[0]-1:Ilats[1],Ilons[0]-1:Ilons[1]]
    #im1  = m.pcolor(lons2, lats2, sst2, cmap='bwr', vmin=cmin,vmax=cmax)
    im1  = m.pcolor(lons2, lats2, cs2, cmap='jet', vmin=cmin,vmax=cmax)
    im2  = m.quiver(lons2[::afreq,::afreq],lats2[::afreq,::afreq],u2[::afreq,::afreq],v2[::afreq,::afreq])
    im3  = polygon_patch(m,ax)
    im4  = ax.quiver(wx,wy,u_w,v_w,scale = 4.0, scale_units='inches',width=.003, color = 'darkgreen') 
    #im5  = ax.Circle((wx,wy),1/3,ec = 'k')
    tx_str = str(year) + ' - ' + mon + ' - ' + str(day).zfill(2)
    tx.set_text(tx_str)
    print i
    print tx_str

if MEEP == 1:
   ani = animation.FuncAnimation(fig, updatefig,frames=np.sum(ndays), blit=False)
  # ani = animation.FuncAnimation(fig, updatefig,frames=8, blit=False)
   ani.save('test_1999.gif', writer = 'imagemagick',fps=4)

else:
   ## MEAN
   u, v,cur_speed = get_vels2()
   u2    = u[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
   v2    = v[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
   cs2   = cur_speed[Ilats[0]-1:Ilats[1],Ilons[0]-1:Ilons[1]]

   #sst = get_sst()
   #sst[:,:] = .35 
   u_w, v_w = get_wind(np.array(range(153,275)))
   #pcol = m.pcolor(lons, lats, sst, cmap='bwr', vmin=cmin,vmax=cmax)
   pcol = m.pcolor(lons2, lats2, cs2, cmap='jet', vmin=cmin,vmax=cmax)
   m.quiver(lons2[::afreq,::afreq],lats2[::afreq,::afreq],u2[::afreq,::afreq],v2[::afreq,::afreq])
   polygon_patch(m,ax)
   m.drawmeridians([120.4,121.4])
   ax.quiver(wx,wy,u_w,v_w,scale = 4.0, scale_units='inches',width=.003,color = 'darkgreen' )
   m.drawparallels([13.4,13.9])
#   ax.quiver(wx,wy,u_w,v_w,headaxislength=3,scale = 4.0, scale_units='inches',width=.004)
#   plt.colorbar(pcol,orientation = 'horizontal')
   #plt.savefig('CIRCLE_1997_cur.png')
   plt.savefig('empty.png')
   #plt.savefig('AVERAGE_1999_cur.png')
