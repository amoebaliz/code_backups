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

def get_dist(tval,temp):
    store_val = np.zeros((temp.shape[1],temp.shape[2]))
    store_val_2 = np.zeros((temp.shape[0]))
    denom = ((sN[tval]-sS[tval])*(sE[tval]-sW[tval]))
    for snt in range(temp.shape[0]):
        temp2 = np.squeeze(temp[snt,:,:])
        avg_T = np.mean(temp2)
        std_T = np.std(temp2)
        val = np.array(np.where(temp2 < (avg_T - 1*std_T)))
        
        #plt.figure()
        #plt.plot(val[1,:],val[0,:],'o')
        #plt.plot([sW[tval],sE[tval]],[sN[tval],sN[tval]],'-k')
        #plt.plot([sE[tval],sE[tval]],[sN[tval],sS[tval]],'-k')
        #plt.plot([sE[tval],sW[tval]],[sS[tval],sS[tval]],'-k')
        #plt.plot([sW[tval],sW[tval]],[sS[tval],sN[tval]],'-k')
        #plt.show()

        for npt in range(val.shape[1]):
            store_val[val[0,npt],val[1,npt]]+=1
            if sN[tval] > val[0,npt] and \
               sS[tval] < val[0,npt] and \
               sE[tval] > val[1,npt] and \
               sW[tval] < val[1,npt]:

               store_val_2[snt]+=1

    per_temp = store_val/(temp.shape[0]) 

    return per_temp, store_val_2

def get_cur(tval,u,v,tndx):
    # interpolation
    u2 = (u[:,:,:-1]+u[:,:,1:])/2
    v2 = (v[:,:-1,:]+v[:,1:,:])/2      

    # composite   
    u3 = np.mean(u2[tndx,:,:],axis=0) 
    v3 = np.mean(v2[tndx,:,:],axis=0)

    # projection
    u4 = u3*(np.cos(angs)) + v3*(np.cos(np.pi/2 + angs))
    v4 = u3*(np.sin(angs)) + v3*(np.sin(np.pi/2 + angs))

    return u4, v4
# ------------------------------------- #
# VIP REGION 
Ilats = [60,232]
Ilons = [168,335]

# SUB-COOLING REGIONS
sN = [100,140,75,120]
sS = [90,100,40,110]
sE = [76,166,70,100]
sW = [47,118,5,80]

# VIP SPECS
VIP = pyroms.grid.get_ROMS_grid('VIP')
lats = VIP.hgrid.lat_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
lons = VIP.hgrid.lon_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
angs = VIP.hgrid.angle_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]

# PLOT SPECS
deg_off = 0.02
cval = .7
afreq = 5
rnud = 0.011
filestr = 'd.png'

fig, ax = plt.subplots(4,2,figsize=(15,16),sharey=True,sharex=True)

yrs = np.arange(1996,1999+1)
for nt in range(len(yrs)):

    ncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/VIP_SST_JJAS_' + str(yrs[nt]) + '.nc'
    uncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/VIP_U_JJAS_' + str(yrs[nt]) + '.nc'
    vncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/VIP_V_JJAS_' + str(yrs[nt]) + '.nc' 
    zncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/VIP_SSH_JJAS_' + str(yrs[nt]) + '.nc'

    fid  = nc.Dataset(ncfile)
    fidu = nc.Dataset(uncfile) 
    fidv = nc.Dataset(vncfile)
    fidz = nc.Dataset(zncfile)

    temp = np.squeeze(fid.variables['temp'][:])
    u    = np.squeeze(fidu.variables['u'][:])
    v    = np.squeeze(fidv.variables['v'][:])
    ssh  = np.squeeze(fidz.variables['zeta'][:]) 

    plt_temp, idx_wts = get_dist(nt,temp)
    tndx = np.where(idx_wts>0)[0] 
    plt_u, plt_v = get_cur(nt,u,v,tndx)
    if nt == 2:
       print tndx
    for ii in range(2):
        m = Basemap(llcrnrlon=np.min(lons)+.05,llcrnrlat=np.min(lats)+.3,urcrnrlon=np.max(lons)+deg_off,urcrnrlat=np.max(lats)-.19,resolution='f',ax =ax[nt,ii] )
        if ii == 0:
           C = m.pcolor(lons, lats, 100*plt_temp,vmin=0,vmax=100*cval,cmap='Blues',zorder=-3)
           # plotting box bounds
           m.plot([lons[sN[nt],sW[nt]],lons[sN[nt],sE[nt]]],[lats[sN[nt],sW[nt]],lats[sN[nt],sE[nt]]],'-w',lw=2)
           m.plot([lons[sN[nt],sE[nt]],lons[sS[nt],sE[nt]]],[lats[sN[nt],sE[nt]],lats[sS[nt],sE[nt]]],'-w',lw=2)
           m.plot([lons[sS[nt],sE[nt]],lons[sS[nt],sW[nt]]],[lats[sS[nt],sE[nt]],lats[sS[nt],sW[nt]]],'-w',lw=2)
           m.plot([lons[sS[nt],sW[nt]],lons[sN[nt],sW[nt]]],[lats[sS[nt],sW[nt]],lats[sN[nt],sW[nt]]],'-w',lw=2)
           # year labels
           ax[nt,ii].text(121.35,13.85, str(yrs[nt]), va = 'top', ha='right', fontsize=16)
        else:
           C2 = m.pcolor(lons, lats, 100*(np.mean(ssh[tndx,:,:],axis=0)-np.min(np.mean(ssh[tndx,:,:],axis=0))),vmin=0,vmax=6, zorder=-3)
           Q  = m.quiver(lons[::afreq,::afreq],lats[::afreq,::afreq],plt_u[::afreq,::afreq],plt_v[::afreq,::afreq])

        ax[nt,ii].plot([lons[0,0],lons[-1,0]],[lats[0,0],lats[-1,0]],'--k',linewidth = 1,zorder=-2)
        ax[nt,ii].plot([lons[0,-1],lons[-1,-1]],[lats[0,-1],lats[-1,-1]],'--k',linewidth = 1,zorder = -2)  
        polygon_patch(m,ax[nt,ii])

        if (nt == 0 and ii == 0):
           m.drawmeridians([120.4,121.4],labels=[False,False,True,False],fontsize=14)
           m.drawparallels([13.4,13.9],labels=[True,False,False,False],fontsize=14)
        else:
           m.drawmeridians([120.4,121.4])
           m.drawparallels([13.4,13.9])
        
        if (nt == 0 and ii == 1):
           qk = ax[nt,ii].quiverkey(Q,.8,.8,.5,r'$0.5  \frac{m}{s}$', labelpos='W',fontproperties={'size':16})

        ax[nt,ii].set_ylim(np.min(lats)+.303,np.max(lats)-.2)
        ax[nt,ii].set_xlim(np.min(lons)+.06,np.max(lons)+deg_off-.02)

    fid.close()
    fidu.close()
    fidv.close()
    fidz.close()
 
plt.subplots_adjust(wspace = 0.05,hspace = 0.05)

cbar_ax  = fig.add_axes([ax[nt,0].get_position().bounds[0]+rnud,0.06,ax[nt,0].get_position().bounds[2]-2*rnud,0.02])
cbar_ax2 = fig.add_axes([ax[nt,1].get_position().bounds[0]+rnud,0.06,ax[nt,1].get_position().bounds[2]-2*rnud,0.02])
plt.colorbar(C,cax=cbar_ax,orientation="horizontal",ticks = [0,35,70])
plt.colorbar(C2,cax=cbar_ax2,orientation="horizontal",ticks = [0,3,6])
cbar_ax.text(.5,.5,'Percent',va='center',ha='center',fontsize=14)
cbar_ax2.text(.5,.5,'Relative SSH (cm)',va='center',ha='center',fontsize=14)
cbar_ax.tick_params(labelsize=14) 
cbar_ax2.tick_params(labelsize=14) 
plt.savefig(filestr)
plt.show()
