import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gsp
from mpl_toolkits.basemap import Basemap
import pyroms as py

mon_ab = ['JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ','DJF']

def make_fig(i,idx):
    fig, ax = plt.subplots(1, figsize = (10,8))
    fig.tight_layout(pad=2)
    plt.subplots_adjust(bottom=0.1, right=0.85, top=0.9, left=0.1)
    tit_str = 'Average ' + mon_ab[i] + ' Camotes Surface Currents'
    plt.title(tit_str, y = 1.04, fontsize=16)
    cbar_ax = fig.add_axes([0.88,0.12,0.015,0.7])
    ax.set_axis_bgcolor((0.7,0.7,0.7))
    ax.set_xlabel('I Direction', labelpad=20, fontsize=14)
    ax.set_ylabel('J Direction', labelpad=30, fontsize=14)

#idx = [0,1,2,12,13,14,24,25,26,36,37,38]
#idx = np.array(idx)+ 4 
#    afreq = 10
    clev = .9
#    levs = np.linspace(0,clev,100,endpoint=True)
    u_store = np.empty([372,252,301])
    v_store = np.empty([372,251,302])
    for nt in range(len(idx)):
        c = idx[nt]+1
        ncfile = 'temp_' + str(c).zfill(5)+'.nc'
        fid = nc.Dataset(ncfile, 'r')
        U = np.squeeze(fid.variables['u'][:])
        V = np.squeeze(fid.variables['v'][:])

        if nt == 0:
           a = 0
        else: 
           a = b
        b = U.shape[0]+a
        print a
        print b
        u_store[a:b,:,:] = U
        v_store[a:b,:,:] = V
    # INTERPOLATE TO RHO-POINTS 
    u_store = u_store[:b,:,:]
    v_store = v_store[:b,:,:]
    print u_store.shape
    u = (u_store[:,1:-1,0:-1]+u_store[:,1:-1,1:])/2
    v = (v_store[:,0:-1,1:-1]+v_store[:,1:,1:-1])/2
    u2 = U

    u = np.float64(u)
    v = np.float64(v)
    # CALCULATE CURRENT SPEED 
    current_speed = np.var(np.sqrt(np.square(u)+np.square(v)),axis=0)
    current_speed = np.ma.masked_where(current_speed>10000,current_speed)
    print np.max(current_speed)
    print np.min(current_speed)
    fid.close()
    x = np.arange(0,current_speed.shape[1])
    y = np.arange(0,current_speed.shape[0])
    X, Y = np.meshgrid(x, y)
    C = ax.pcolor(X,Y,current_speed,vmin = 0, vmax=.9)
    fig.colorbar(C,cax=cbar_ax,ticks = [0, clev/4, clev/2, 3*clev/4, clev])
    cbar_ax.text(-1.2,1.06,'Current Speed \n Variance \n (m/s)',multialignment='center')

    u_new = np.mean(u,axis=0)
    v_new = np.mean(v,axis=0)

#    Q = ax.quiver(X[::afreq,::afreq],Y[::afreq,::afreq],u_new[::afreq,::afreq],v_new[::afreq,::afreq],width = 0.0015)
#    qk = ax.quiverkey(Q,.1,.94,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})
    save_tit = 'camotes_var_'+mon_ab[i]+'.png'
    plt.savefig(save_tit)
# ITERATE OVER ALL MONTH COMBOS

idx_vals = np.empty([4,14])
idx_vals[:] = np.NAN
idx_vals2 = np.array(range(48)).reshape((4,12))

idx_vals[:,:10] = idx_vals2[:,:10]
idx_vals[:-1,10:12]=idx_vals2[:-1,10:]
idx_vals[:,12:14]= idx_vals2[:,:2]

for nt in range(12):
    idx = np.squeeze(np.ravel(idx_vals[:,nt:nt+3]))
    idx = idx[~np.isnan(idx)]
    idx2 = [int(idx_val) for idx_val in idx]
    print idx2
    make_fig(nt,idx2)

