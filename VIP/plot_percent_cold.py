import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyroms
from matplotlib.collections import PolyCollection

def polygon_patch(mapid,axs):
#    mapid.drawmapboundary(fill_color=[.9,.97,1])
    mapid.drawcoastlines(linewidth=0)
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black', \
         facecolor='0.8', closed=False)
    axs.add_collection(lc)

def get_dist(temp):
    store_val = np.zeros((temp.shape[1],temp.shape[2]))
    for nt in range(temp.shape[0]):
        temp2 = np.squeeze(temp[nt,:,:])
        avg_T = np.mean(temp2)
        std_T = np.std(temp2)
        val = np.array(np.where(temp2 < (avg_T - 1*std_T)))
        for npt in range(val.shape[1]):
            store_val[val[0,npt],val[1,npt]]+=1
    per_temp = store_val/(temp.shape[0])
 
    return per_temp
# -------------------------------------#
VIP = pyroms.grid.get_ROMS_grid('VIP')
Ilats = [60,232]
Ilons = [168,335]
deg_off = 0.02
cval = .7
yrs = np.arange(1996,1999+1)
fig, ax = plt.subplots(2,2,figsize=(15,8),sharey=True,sharex=True)
cbar_ax = fig.add_axes([0.16,0.04,0.7,0.04])
lats = VIP.hgrid.lat_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
lons = VIP.hgrid.lon_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1] 
for nt in range(len(yrs)): 
    i = (nt+1)/2 - nt/2
    j = nt/2
#    ncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(yrs[nt]) + '/VIP_normSST_MJJAS_' + str(yrs[nt]) + '.nc'
    ncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/VIP_SST_JJAS_' + str(yrs[nt]) + '.nc'
    fid = nc.Dataset(ncfile)
    temp = np.squeeze(fid.variables['temp'][:])
    plt_temp = get_dist(temp)
#    plt_temp = np.mean(temp,axis = 0)
#    C = ax[nt].pcolor(plt_temp,vmin=-1*cval,vmax=cval,cmap='bwr')
    m = Basemap(llcrnrlon=np.min(lons)+.05,llcrnrlat=np.min(lats)+.3,urcrnrlon=np.max(lons)+deg_off,urcrnrlat=np.max(lats)-.19,resolution='f',ax =ax[j,i] )
    C = m.pcolor(lons, lats, 100*plt_temp,vmin=0,vmax=100*cval,cmap='Blues',zorder=-3)

    ax[j,i].plot([lons[0,0],lons[-1,0]],[lats[0,0],lats[-1,0]],'--k',linewidth = 0.5,zorder=-2)
    ax[j,i].plot([lons[0,-1],lons[-1,-1]],[lats[0,-1],lats[-1,-1]],'--k',linewidth = 1,zorder = -2)

    polygon_patch(m,ax[j,i])
    ax[j,i].text(121.35,13.85, str(yrs[nt]), va = 'top', ha='right', fontsize=16)
    if nt == 0: 
       m.drawmeridians([120.4,121.4],labels=[False,False,True,False],fontsize=14)
       m.drawparallels([13.4,13.9],labels=[True,False,False,False],fontsize=14)
    else:
       m.drawmeridians([120.4,121.4])
       m.drawparallels([13.4,13.9])
#    ax[j,i].figure.canvas.draw() 
    ax[j,i].set_ylim(np.min(lats)+.303,np.max(lats)-.19)
    ax[j,i].set_xlim(np.min(lons)+.06,np.max(lons)+deg_off-.01)
#    ax[j,i].set_title(str(yrs[nt]))
    fid.close()
plt.subplots_adjust(wspace = 0.05,hspace = 0.05)
plt.colorbar(C,cax=cbar_ax,orientation="horizontal",ticks = [0,35,70])
plt.text(.5,.5,'Percent',va='center',ha='center',fontsize=14)
cbar_ax.tick_params(labelsize=14) 
#fig.tight_layout()
plt.savefig('c.png')
plt.show()
