# COMPARING AVERAGE SS CONDITIONS FOR CTROMS AND VIP
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

ctroms_dir = '/t3/workdir/liz/VIP/VIP_grid/VIP_interpolations/ctroms_vel/CTROMS_surf_2007.nc'
vip_file = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150313/VIP_surf_2007.nc' 
fid = [nc.Dataset(vip_file,'r'), nc.Dataset(ctroms_dir,'r')]

vars = ['temp','salt','u']
cmin = [27,30,0]
cmax = [31,35,.34]
afreq = 10

fig,axarr = plt.subplots(3,2, figsize = (14,14),sharey='row',sharex='col')
plt.tight_layout(h_pad=5)
plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9, left=0.1)
for nv in range(len(vars)):
    cbar_ax = fig.add_axes([0.93,.7-(0.29*(nv)),0.015,0.15]) 
    for ncol in range(0,2):
       var = fid[ncol].variables[vars[nv]][:]
       if nv == 2:
           V = fid[ncol].variables['v'][:]   
           u = np.squeeze((var[:,:,1:-1,0:-1]+var[:,:,1:-1,1:])/2)
           v = np.squeeze((V[:,:,0:-1,1:-1]+V[:,:,1:,1:-1])/2)
           var = np.sqrt(np.square(u)+np.square(v))
           X, Y = np.meshgrid(np.arange(0,var.shape[1]), np.arange(0,var.shape[0]))
           Q = axarr[nv,ncol].quiver(X[::afreq,::afreq],Y[::afreq,::afreq],u[::afreq,::afreq],v[::afreq,::afreq],width = 0.0015,zorder=2)
       axarr[nv,ncol].set_axis_bgcolor((0.7,0.7,0.7))
       C = axarr[nv,ncol].contourf(np.squeeze(var),np.linspace(cmin[nv],cmax[nv],100, endpoint=True),zorder=1)
       axarr[nv,ncol].set_xlim(0,var.shape[1])
       axarr[nv,ncol].yaxis.tick_left()
       axarr[nv,ncol].xaxis.tick_bottom()  
    fig.colorbar(C,cax=cbar_ax,ticks = [cmin[nv],float((cmin[nv]+cmax[nv]))/2, cmax[nv]])
# plot quiver arros
# fix colorbars
# titles
#  


plt.show()


