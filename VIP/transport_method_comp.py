# Calculate the heat content above a particular depth level
import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as pltd

def get_grid_info(grid):
   
    GRD = pyroms.grid.get_ROMS_grid(grid)
 
    z_w = GRD.vgrid.z_w[:][:,tsct_eta,tsct_xi] 
    dx = GRD.hgrid.dx[tsct_eta,tsct_xi]
    dy = GRD.hgrid.dy[tsct_eta,tsct_xi]

    return z_w, dx, dy
#### -----------------

trm_strt = 31 
trm_end =  62
tsct_xi = 278 
tsct_eta = np.array(np.arange(114+trm_strt,114+trm_end))

z_w, dx, dy =  get_grid_info('VIP')

n_dep = -2
ref = dt.datetime(1900,1,1,0,0)

fid = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/u_all_278_96_99.nc')
u_vel = np.ma.mean(fid.variables['u'][:,:,trm_strt:trm_end,:],axis=3)

fid2 = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/SSH_1996-1999.nc')
var = np.squeeze(fid2.variables['zeta'][:])
var = (np.squeeze(np.mean(var[:,65:149,180],axis=1))-np.squeeze(np.mean(var[:,149:199,318],axis=1)))

#fid2 = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/SSH_1996-1999.nc')
#var = np.squeeze(fid2.variables['zeta'][:])
#var = (np.squeeze(np.mean(var[:,65:149,180],axis=1))-np.squeeze(np.mean(var[:,149:199,318],axis=1)))

dep_vals = -1*np.array(np.arange(0,500,1))
cor_vals = np.zeros((len(dep_vals)))
print cor_vals.shape
for nt in range(len(dep_vals)):

    n_dep = dep_vals[nt]
    z_w_1 = np.ma.where(z_w>n_dep,z_w,n_dep)
    dz = z_w_1[1:,:]-z_w_1[:-1,:]
    vol = np.tile(dz*dy,(u_vel.shape[0],1,1))
    transp = vol*u_vel
    cor_vals[nt] = np.corrcoef(np.sum(np.sum(transp,axis=2),axis=1),var)[0,1]  
#pos_t = np.sum(np.sum(np.ma.where(transp>0.0,transp,0),axis=2),axis=1)
#neg_t = np.sum(np.sum(np.ma.where(transp<0.0,transp,0),axis=2),axis=1)

#plt.figure()
#plt.plot(pos_t)
#plt.plot(neg_t)

#fid2 = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/vip_transports/vip_u_2m_transport_col278_rows144-178.nc')
#pos_t2 = fid2.variables['pos_transport'][:]
#neg_t2 =  fid2.variables['neg_transport'][:]

#plt.figure()
#plt.plot(pos_t2)
#plt.plot(neg_t2)

#plt.figure()
#plt.plot(pos_t/np.max(pos_t))
#plt.plot(pos_t2/np.max(pos_t2))

#plt.figure()
#plt.plot(pos_t-pos_t2)
#print np.corrcoef(neg_t,neg_t2)[0,1]
#print np.mean(pos_t/pos_t2)
#print np.min(pos_t/pos_t2)
#print np.max(pos_t/pos_t2)
plt.figure()
plt.plot(dep_vals,cor_vals)
plt.show()

