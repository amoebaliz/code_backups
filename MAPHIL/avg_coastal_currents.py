import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0)
    mapid.drawmapboundary()
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor=(cv,cv,cv), closed=False)
    axs.add_collection(lc)


# FYI - THESE ARE THE LIMITS USED FOR SUB-DIVIDING 
ilims = [278, 300]
jlims = [ 74, 178]
# ncrcat -d s_rho,49 -d eta_u,74,179 -d xi_u,278,300 -d eta_v,74,178 -d xi_v,278,301 -v u,v /t1/scratch/liz/tmpdir_MaPhil-LD.HCo05T/outputs/2010/MaPhil-LD.HCo05T_avg_2010-0[12345]* 1_2010.nc
cv = 0.8
grd = pyroms.grid.get_ROMS_grid('MaPhil')
#lats = grd.hgrid.lat_rho[jlims[0]:jlims[1]+1,ilims[0]:ilims[1]+1]
#lons = grd.hgrid.lon_rho[jlims[0]:jlims[1]+1,ilims[0]:ilims[1]+1]
lats = grd.hgrid.lat_v[jlims[0]:jlims[1]+1,ilims[0]:ilims[1]+2]
lons = grd.hgrid.lon_v[jlims[0]:jlims[1]+1,ilims[0]:ilims[1]+2]
ncfile = '/t3/workdir/liz/MODELS/MAPHIL/PostProc/all_cur.nc'
fid = nc.Dataset(ncfile)

u = np.squeeze(fid.variables['u'][:])
v = np.squeeze(fid.variables['v'][:])
print u.shape
uavg = np.mean(u,axis=0)
vavg = np.mean(v,axis=0)

print vavg.shape
print np.mean(v[:,:,15:])
print np.std(v[:,:,15:])
print np.max(v[:,:,15:])
print np.min(v[:,:,15:])

test_vals = np.mean(np.mean(v,axis=2),axis=1)
#plt.figure()
#plt.plot(test_vals)
#plt.show()

u_interp = (uavg[:-1,:]+uavg[1:,:])/2
v_interp = (vavg[:,:-1]+vavg[:,1:])/2
sped = np.sqrt(np.square(u_interp)+np.square(v_interp))

moff = 0.05
afreq = 5
u_interp[:]=0

fig = plt.figure(figsize=(8,10))
ax = fig.add_subplot(111)
m = Basemap(llcrnrlon=np.min(lons)-moff,llcrnrlat=np.min(lats),urcrnrlon=np.max(lons),urcrnrlat=np.max(lats),resolution='f',ax =ax )
col_val = 0.1
v = np.linspace(-1*col_val,col_val, 50, endpoint=True)
#cs = m.contourf(lons,lats,v_interp,v,iedgecolors='face',cmap='bwr')
cs = m.pcolor(lons,lats,vavg,vmin=-1*col_val,vmax = col_val,cmap='bwr')
m.drawmeridians([124.6,124.8],labels=[0,0,0,1], fontsize=14)
m.drawparallels([10.7,11],labels=[1,0,0,0], fontsize=14)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cs, cax=cax,ticks=[-1*col_val,0,col_val])
cax.set_title('Along-shore Velocity (m/s)', rotation=270, va='center', ha='center', fontsize=16, y=.5)
m.plot([lons[-1,-1],lons[0,-1]],[lats[-1,-1],lats[0,-1]],'-k')

#q = m.quiver(lons[::afreq,::afreq],lats[::afreq,::afreq],u_interp[::afreq,::afreq],v_interp[::afreq,::afreq])
#ax.quiverkey(q,.1,.94,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})
#plt.colorbar(cs)
polygon_patch(m,ax)
plt.show()
