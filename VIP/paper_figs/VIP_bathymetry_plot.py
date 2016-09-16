import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False)
#         facecolor='.8', closed=False)
    axs.add_collection(lc)

def outline_mask(mapid,mask_img,val,x0,y0,x1,y1):
    mapimg = (mask_img == val)
    ver_seg = np.where(mapimg[:,1:] != mapimg[:,:-1])
    hor_seg = np.where(mapimg[1:,:] != mapimg[:-1,:])
  
    l = []
    v = []
    # horizonal segments
    for p in zip(*hor_seg):
        if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
           v.append((lons[p[0]+1,p[1]],lats[p[0]+1,p[1]]))
        else :
           l.append((lons[p[0]+1,p[1]],lats[p[0]+1,p[1]])) 

        if p[1] == mask_img.shape[1] - 1 :
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
               v.append((lons[p[0]+1,p[1]],lats[p[0]+1,p[1]]))         
           else:
               l.append((lons[p[0]+1,p[1]],lats[p[0]+1,p[1]])) 
        else :
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((lons[p[0]+1,p[1]+1],lats[p[0]+1,p[1]+1]))
           else:
              l.append((lons[p[0]+1,p[1]+1],lats[p[0]+1,p[1]+1]))

        l.append((np.nan,np.nan))
        v.append((np.nan,np.nan)) 
    #vertical segments
    for p in zip(*ver_seg):
        if p[1] == mask_img.shape[1]-1:
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((lons[p[0],p[1]],lats[p[0],p[1]]))
              v.append((lons[p[0]+1,p[1]],lats[p[0]+1,p[1]]))
           else:
              l.append((lons[p[0],p[1]],lats[p[0],p[1]]))
              l.append((lons[p[0]+1,p[1]],lats[p[0]+1,p[1]]))    
        else:
           if (vip_eta[0] < p[0] < vip_eta[1] and vip_xi[0] < p[1] < vip_xi[1]):
              v.append((lons[p[0],p[1]+1],lats[p[0],p[1]+1]))
              v.append((lons[p[0]+1,p[1]+1],lats[p[0]+1,p[1]+1])) 
           else:
              l.append((lons[p[0],p[1]+1],lats[p[0],p[1]+1]))
              l.append((lons[p[0]+1,p[1]+1],lats[p[0]+1,p[1]+1]))

        l.append((np.nan, np.nan))
        v.append((np.nan, np.nan))
    segments = np.array(l)
    vip_segments = np.array(v)
    print vip_segments.shape
    mapid.plot(segments[:,0], segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+2) 
#    mapid.plot(vip_segments[:,0], vip_segments[:,1], latlon=True, color=(0,0,1), linewidth=.75,zorder=map_order+3)
    mapid.plot(vip_segments[:,0], vip_segments[:,1], latlon=True, color=(0,0,0), linewidth=.75,zorder=map_order+3)
def plot_transects(mapid,z_order):
    for ntrans in range(len(tsect_xi)):
        eta = tsect_eta[2*ntrans:(2*ntrans+2)]
        xi  = tsect_xi[ntrans]

        t_lats = lats[eta,(xi, xi)]
        t_lons = lons[eta,(xi, xi)]

        tid_lats = lats[eta[0]+eta_offsets[ntrans],xi+xi_offsets[ntrans]] 
        tid_lons = lons[eta[0]+eta_offsets[ntrans],xi+xi_offsets[ntrans]]
 
        mapid.plot(t_lons, t_lats, color='w', lw=1, ls='--',dashes=(3,2), zorder =z_order)
        mapid.plot(tid_lons, tid_lats, 'o', ms=20, mfc='w', mew=1, alpha=0.9, zorder =z_order)
        plt.text(tid_lons, tid_lats, t_id[ntrans],fontsize=14,ha = 'center',va ='center',zorder =z_order+1)

def plot_stations(mapid,z_order):
    for nt in range(len(vip_stations)/2):
        st_lat = lats[vip_stations[2*nt+1],vip_stations[2*nt]]
        st_lon = lons[vip_stations[2*nt+1],vip_stations[2*nt]]
        mapid.plot(st_lon,st_lat,'D',ms=8, mfc=(1,1,0),mec='k',mew=1, zorder=z_order)

# ---------------------------------------------------- #

VIP      = pyroms.grid.get_ROMS_grid('VIP')
bathy    = VIP.vgrid.h
mask_rho = VIP.hgrid.mask_rho
lats     = VIP.hgrid.lat_rho
lons     = VIP.hgrid.lon_rho

mask_val = -50
bathy[mask_rho==0] = mask_val
mask_plot = np.ma.masked_where(mask_rho>0,mask_rho)
vip_eta = [60,240]
vip_xi  = [177,338]
# Transect Info
tsect_eta = np.array([65,149,102,152,145,178,146,199])
tsect_xi = np.array([180,240,280,318])
xi_offsets = [8,5,3,7]
eta_offsets = [-8,-12,-18,-11]
t_id = ['i','ii','iii','iv']


# Station Info
vip_stations = [126, 320] # 281, 156]  
# (Manila, Transect iii)

# FIGURE DETAILS
fig = plt.figure(figsize=(18,12))
ax1 = fig.add_subplot(111)

# Figure Domain
m_offset = 0.1
lats1 = [lats[0,-1] -m_offset, lats[-1, 0] + m_offset]
lons1 = [lons[0, 0] -m_offset, lons[-1,-1] + m_offset]
m = Basemap(llcrnrlon=lons1[0],llcrnrlat=lats1[0],urcrnrlon=lons1[1],urcrnrlat=lats1[1],resolution='f')

map_order = 30
P = m.pcolormesh(lons,lats,bathy,vmin=0,vmax=3600,edgecolors='face',cmap='nipy_spectral_r',zorder=map_order)
P.cmap.set_under('w')

# MAP DETAILING
outline_mask(m,bathy,mask_val,lons[0,0],lats[0,0],lons[-1,-1],lats[-1,-1])
plot_transects(m,map_order+3)
plot_stations(m,map_order+3)


# DOMAIN OUTLINE
m.plot((lons[0,0],lons[-1,0]),(lats[0,0],lats[-1,0]),linewidth=2,color='k',zorder=map_order+1)
m.plot((lons[-1,0],lons[-1,-1]),(lats[-1,0],lats[-1,-1]),linewidth=2,color='k',zorder=map_order+1)
m.plot((lons[-1,-1],lons[0,-1]),(lats[-1,-1],lats[0,-1]),linewidth=2,color='k',zorder=map_order+1)
m.plot((lons[0,-1],lons[0,0]),(lats[0,-1],lats[0,0]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax1)
m.drawmeridians([120,122], labels=[0,0,1,0], fmt='%d', fontsize=18)
m.drawparallels([13,15], labels=[1,0,0,0], fmt='%d', fontsize=18)

ax1.figure.canvas.draw()
pos1 = ax1.get_position()
cbar_ax = fig.add_axes([pos1.x1+.01,pos1.y0,0.025,pos1.y1-pos1.y0])
cbar = plt.colorbar(P,cax=cbar_ax,ticks=[0,1200,2400,3600])
cbar.ax.invert_yaxis() 
cbar_ax.tick_params(labelsize=14) 
cbar_ax.set_title('Depth (m)', y=0.5, rotation=270, fontsize = 18, va='center',ha = 'center')

ax1.text(121.2,14.08,'Luzon',fontsize = 18, ha='center', va='center',zorder = map_order+3)
ax1.text(121.08,13.1,'Mindoro',fontsize = 18, ha='center', va='center',zorder = map_order+3)
#ax1.text(121.0,14.53,'Manila',style='italic', va='center',zorder = map_order+3)

ax1.text(121.982,13.402,'Marinduque', fontsize = 11, ha = 'center', va='center',zorder = map_order+3)
ax1.text(120.17,13.8,'Lubang', fontsize = 8, ha = 'center', va='center',zorder = map_order+3)

ax1.text(121.548,13.42,'Tayabas \n Bay',fontsize = 18, color = 'w', ha='center', va='center',zorder = map_order+3)
ax1.text(119.71,13.2,'South China \n Sea',fontsize = 18, ha='center', va='center',zorder = map_order+3)

ax1.text(120.8,13.77,'Balayan \n Bay', color = 'w', fontsize = 11, ha = 'center', va='center',zorder = map_order+3)
ax1.text(120.98,13.68,'*', color = 'w', fontsize = 16, ha = 'center', va='center',zorder = map_order+3)

#plt.title('Verde Island Passage\n Model Domain and Bathymetry (meters)',fontsize=18)
#plt.tight_layout()
plt.savefig('blank.png')
plt.show()




