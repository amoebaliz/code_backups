import numpy as np
import netCDF4 as nc
import datetime as dt
import pyroms
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs,map_order):
    mapid.drawcoastlines(linewidth=0)
    mapid.drawmapboundary(fill_color=[.9,.97,1])
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(1,1,1), closed=False,zorder=map_order+3)
    axs.add_collection(lc)

##############################################
VIP      = pyroms.grid.get_ROMS_grid('VIP')
bathy    = VIP.vgrid.h
mask_rho = VIP.hgrid.mask_rho
lats     = VIP.hgrid.lat_rho
lons     = VIP.hgrid.lon_rho

dhw_file = '/t3/workdir/liz/MODELS/VIP/PostProc/DHW/Max_VIP_DHW_CoRTclim_1998.nc'
fid = nc.Dataset(dhw_file)

cmap_file = '/home/frederic/python/cmap/dhw_noaa.cmap'
DHW = fid.variables['max_DHW'][:]
DHW = np.ma.masked_where(DHW<0,DHW)
cmap = np.loadtxt(cmap_file)
lev_seven = np.array((159,212,0))
cmap2 = np.vstack((cmap[0:7,:],lev_seven,cmap[7:,:]))
cmap = cmap2/256.
dhw_noaa = colors.ListedColormap(cmap)

mask_val = -50
bathy[mask_rho==0] = mask_val
mask_plot = np.ma.masked_where(mask_rho>0,mask_rho)
vip_eta = [60,240]
vip_xi  = [177,338]

# FIGURE DETAILS
fig = plt.figure(figsize=(15,12))
ax1 = fig.add_subplot(111)

# Figure Domain
m_offset = 0.1
lats1 = [lats[0,-1] -m_offset, lats[-1, 0] + m_offset]
lons1 = [lons[0, 0] -m_offset, lons[-1,-1] + m_offset]
m = Basemap(llcrnrlon=lons1[0],llcrnrlat=lats1[0],urcrnrlon=lons1[1],urcrnrlat=lats1[1],resolution='f')

map_order = 30
cmin=0
cmax=16
P = m.pcolormesh(lons,lats,np.squeeze(DHW),vmin=cmin,vmax=cmax,edgecolors='face',cmap=dhw_noaa,zorder=map_order)
P.cmap.set_under('w')

# DOMAIN OUTLINE
m.plot((lons[0,0],lons[-1,0]),(lats[0,0],lats[-1,0]),linewidth=2,color='k',zorder=map_order+4)
m.plot((lons[-1,0],lons[-1,-1]),(lats[-1,0],lats[-1,-1]),linewidth=2,color='k',zorder=map_order+4)
m.plot((lons[-1,-1],lons[0,-1]),(lats[-1,-1],lats[0,-1]),linewidth=2,color='k',zorder=map_order+4)
m.plot((lons[0,-1],lons[0,0]),(lats[0,-1],lats[0,0]),linewidth=2,color='k',zorder=map_order+4)

polygon_patch(m,ax1,map_order)
m.drawmeridians([120,122], labels=[0,0,0,1], fmt='%d', fontsize=18)
m.drawparallels([13,15], labels=[1,0,0,0], fmt='%d', fontsize=18)
cbar = plt.colorbar(P)
#cbar.ax.invert_yaxis()
plt.title('Verde Island Passage\n Max 1998 DHW',fontsize=18)

plt.show()
