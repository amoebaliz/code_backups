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
#         facecolor=(1,1,1), closed=False)
          facecolor='.8', closed=False)
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
#        mapid.plot(tid_lons, tid_lats, 'o', ms=19, mfc='w', mew=1, alpha=0.9, zorder =z_order)
#        plt.text(tid_lons, tid_lats, t_id[ntrans],fontsize=13,ha = 'center',va ='center',zorder =z_order+1)

def plot_stations(mapid,z_order):
 
    for nt in range(len(vip_stations)/2):
#        st_lat = lats[vip_stations[2*nt+1],vip_stations[2*nt]]
#        st_lon = lons[vip_stations[2*nt+1],vip_stations[2*nt]]
        mapid.plot(vip_stations[2*nt],vip_stations[2*nt+1],'o',ms=6, mfc=(1,1,0),mec='k',mew=1, zorder=z_order)

# ---------------------------------------------------- #

VIP      = pyroms.grid.get_ROMS_grid('MaPhil')
bathy    = VIP.vgrid.h
mask_rho = VIP.hgrid.mask_rho
lats     = VIP.hgrid.lat_rho
lons     = VIP.hgrid.lon_rho

mask_val = -50
bathy[mask_rho==0] = mask_val
mask_plot = np.ma.masked_where(mask_rho>0,mask_rho)
#vip_eta = [60,240]
#vip_xi  = [177,338]


# Station Info
vip_stations = [124.6293029, 10.91280546, \
124.6608268, 10.88025355, \
124.6820712, 10.83845005, \
124.7142805, 10.80315692, \
124.7372381, 10.76409463, \
124.7478603, 10.71920727, \
124.7468324, 10.67294929, \
124.7266159, 10.63183109, \
124.7197629, 10.59002759, \
124.7125672, 10.54376961, \
124.6200513, 10.8076114, \
124.6358132, 10.67603316, \
124.5693388, 10.56364341, \
124.786665,  10.743685]

# (Manila, Transect iii)

# FIGURE DETAILS
fig = plt.figure(figsize=(15,12))
ax1 = fig.add_subplot(111)

# Figure Domain
m_offset = 0.1
#lats1 = [lats[0,-1] -m_offset, lats[-1, 0] + m_offset]
lons1 = [lons[0, 0] -m_offset, lons[-1,-1] + m_offset]
m = Basemap(llcrnrlon=np.min(lons),llcrnrlat=np.min(lats),urcrnrlon=np.max(lons),urcrnrlat=np.max(lats),resolution='f')

map_order = 0
P = m.pcolormesh(lons,lats,bathy,vmin=0,vmax=3600,edgecolors='face',cmap='nipy_spectral_r',zorder=map_order)
P.cmap.set_under('w')

# MAP DETAILING
# outline_domain(m)
#outline_mask(m,bathy,mask_val,lons[0,0],lats[0,0],lons[-1,-1],lats[-1,-1])
#plot_transects(m,map_order+3)
plot_stations(m,map_order+3)


# DOMAIN OUTLINE
m.plot((lons[0,0],lons[-1,0]),(lats[0,0],lats[-1,0]),linewidth=2,color='k',zorder=map_order+1)
m.plot((lons[-1,0],lons[-1,-1]),(lats[-1,0],lats[-1,-1]),linewidth=2,color='k',zorder=map_order+1)
m.plot((lons[-1,-1],lons[0,-1]),(lats[-1,-1],lats[0,-1]),linewidth=2,color='k',zorder=map_order+1)
m.plot((lons[0,-1],lons[0,0]),(lats[0,-1],lats[0,0]),linewidth=2,color='k',zorder=map_order+1)

polygon_patch(m,ax1)
#m.drawmeridians([120,122], labels=[0,0,0,1], fmt='%d', fontsize=18)
#m.drawparallels([13,15], labels=[1,0,0,0], fmt='%d', fontsize=18)
#cbar = plt.colorbar(P, orientation='horizontal')
#cbar.ax.invert_yaxis() 
plt.title('Camotes Sea\n Model Domain and Bathymetry (meters)',fontsize=18)
plt.tight_layout()
plt.savefig('blank.png')
plt.show()




