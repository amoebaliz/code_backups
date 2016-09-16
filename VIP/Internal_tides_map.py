import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
import matplotlib.gridspec as gsp
from matplotlib.ticker import MultipleLocator

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0)
    mapid.drawmapboundary()
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(cv,cv,cv), closed=False)
    axs.add_collection(lc)

def draw_transect_lines(ax):

    for nt in range(len(istart)):
        i0 = grd.hgrid.lon_rho[jstart[nt],istart[nt]] 
        iN = grd.hgrid.lon_rho[jend[nt],iend[nt]]
        j0 = grd.hgrid.lat_rho[jstart[nt],istart[nt]]
        jN = grd.hgrid.lat_rho[jend[nt],iend[nt]]
        
        m.plot([i0,iN],[j0,jN],'-k')

        # Plot orientation points
        if nt == 0:
           px = i0; py = j0
        else:
           px = iN; py = jN
        m.plot(px,py,'o',ms=4,mfc='k')

def get_dist(i):

    dist=np.sqrt(np.square((np.sum(dx[jstart[i],istart[i]:iend[i]+1])   + \
                            np.sum(dx[jend[i],istart[i]:iend[i]+1]))/2) + \
                 np.square((np.sum(dy[jstart[i]:jend[i]+1,istart[i]])   + \
                            np.sum(dy[jstart[i]:jend[i]+1,iend[i]]))/2))
    return dist

def plot_map(ax):
    global m, v, cs
    m = Basemap(llcrnrlon=mer[0]-moff-.1,llcrnrlat=par[0]-moff,urcrnrlon=mer[1]+moff+.1,urcrnrlat=par[1]+moff,resolution='h',ax =ax )
    map_order = 0
    v = np.linspace(-1*col_val,col_val, 50, endpoint=True)

    cs = m.contourf(lons,lats,1000*w_dept,v,iedgecolors='face',extend = 'both', cmap='bwr')
    cs.cmap.set_over('r')
    cs.cmap.set_under('b')
    polygon_patch(m,ax)
    m.drawmeridians(mer,labels=[0,0,0,1],fontsize=14)
    m.drawparallels(par,labels=[1,0,0,0], fontsize=14)
    #for nt in range(len(par)):
        
       # txp = plt.text(np.mean(mer),par[nt],str(par[nt])+ u'\N{DEGREE SIGN}N',fontsize=16,ha='center',va='center',zorder=map_order+3)
       # txp.set_bbox(dict(color=[cv,cv,cv], alpha=0.7))    
       # txm = plt.text(mer[nt],np.mean(par),str(mer[nt])+ u'\N{DEGREE SIGN}E',fontsize=16,ha='center',va='center',zorder=map_order+3)
       # txm.set_bbox(dict(color='w', alpha=0.7))
    draw_transect_lines(ax)
def plot_transect(ax,i):

        transect, z, lon, lat = pyroms.tools.transect(var[:-1,:,:], istart[i], iend[i], jstart[i], jend[i], grd,Cpos='rho')
        x = np.tile(np.linspace(0, 1, num=transect[:].shape[1], endpoint=True),(transect[:].shape[0],1))
        ax.set_axis_bgcolor((cv, cv, cv))
        cs1 = ax.contourf(x,-1*z,1000*transect,v,extend = 'both', cmap='bwr')
 
        cs1.cmap.set_over('r')
        cs1.cmap.set_under('b')
        ax.set_ylim(0,zlim[nl-1])
        ax.set_yticks([0,tmx[nl-1]/2, tmx[nl-1]])
        ax.xaxis.set_minor_locator(minorLocator)
        ax.tick_params(axis='both', which='both', top='off',labelsize=14)
        ax.plot([0,1],[zlim[nl-1],zlim[nl-1]],'k',lw = 2)
        ax.invert_yaxis()

        dist=np.sqrt(np.square((np.sum(dx[jstart[i],istart[i]:iend[i]+1])   + \
                                np.sum(dx[jend[i],istart[i]:iend[i]+1]))/2) + \
                     np.square((np.sum(dy[jstart[i]:jend[i]+1,istart[i]])   + \
                                np.sum(dy[jstart[i]:jend[i]+1,iend[i]]))/2))
        print dist

# ---------------------------------------------------- #
grd = pyroms.grid.get_ROMS_grid('VIP')
dx  = grd.hgrid.dx
dy  = grd.hgrid.dy

Ilats = [10,600]
Ilons = [110,460]

lats = grd.hgrid.lat_rho
lons = grd.hgrid.lon_rho

# TRANSECTS
istart = np.array([220,320])
iend   = np.array([260,400])
jstart = np.array([105,190])
jend   = np.array([133,252])

depth_val = 50

ncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/1998/VIP-LD.HCo13T_avg_1998-07-24T00:00:00.nc'
fid = nc.Dataset(ncfile)
var = np.squeeze(fid.variables['w'][:])
w_dept,x,y = pyroms.tools.zslice(var, depth_val, grd, Cpos='w', vert=False, mode='linear')

col_val = 1
moff = 0.15

# FIGURE DETAILS
nrows=4
cv = 0.8
mer = [120.5,121.5]
par = [13.4,13.9]
label = ['a','b','c']
zlim = [800,650]
tmx = [800,600] 
minorLocator = MultipleLocator(.25)
fig = plt.figure(figsize=(8,10))
gs = gsp.GridSpec(nrows,1)
gs.update(bottom = 0.1, left = 0.15, right = 0.82, top = 0.9, hspace=.3)

# MAP
ax = fig.add_subplot(gs[0:nrows-2])
plot_map(ax)
ax.figure.canvas.draw()
pos1 = ax.get_position()

# TRANSECT PLOTS
ax1 = fig.add_subplot(gs[nrows-2])
ax2 = fig.add_subplot(gs[nrows-1],sharex=ax1)
nl = 0
for ax in fig.get_axes():           
    if nl > 0:
       plot_transect(ax,nl-1)
       ax.figure.canvas.draw()
       pos = ax.get_position() 
    else:
       pos=pos1

    ax.text(0.02, .95, label[nl], transform=ax.transAxes,fontsize=16,va='top')
    nl+=1

ax1.tick_params(labelbottom='off')
ax1.plot(0,zlim[0],'o',ms=8,mfc='k',mew=1,clip_on=False,zorder=100)
ax2.plot(1,zlim[1],'o',ms=8,mfc='k',mew=1,clip_on=False,zorder=100)
ax2.set_xticks([0,0.5,1.0])
ax2.set_xlabel('Fraction of Transect Length',fontsize=16,labelpad=30)
ax2.set_ylabel('Depth (m)',y=1.1,fontsize=16,labelpad=30)
pos2 = ax2.get_position()
cbar_ax = fig.add_axes([pos1.x1+.05,pos2.y0,0.05,pos1.y1-pos2.y0])
cbar_ax.tick_params(labelsize=14)
cbar_ax.set_title('Vertical Velocity (mm/s)', rotation=270, va='center', ha='center', fontsize=16, y=.5)

plt.colorbar(cs,cax=cbar_ax,ticks=[-1,0,1])
plt.show()
