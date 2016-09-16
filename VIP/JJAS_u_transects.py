import netCDF4 as nc
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime as dt
import pyroms

def get_vels(ncfile):
     fid = nc.Dataset(ncfile,'r')
     u = np.squeeze(fid.variables['u'][:,:,130:200,278],axis=0)
     u = np.ma.masked_where(u>1000,u)
     zeta = np.squeeze(fid.variables['zeta'][:])
     z = pyroms.vgrid.z_w(VIP.vgrid.h, VIP.vgrid.hc, VIP.vgrid.N, VIP.vgrid.s_rho, VIP.vgrid.Cs_r, zeta, VIP.vgrid.Vtrans)
     z_vals = np.squeeze(z[:,:,130:200,290])
     z_vals = np.ma.masked_where(z_vals>0,z_vals)

     return u, z_vals

VIP = pyroms.grid.get_ROMS_grid('VIP')
mon_days = [31,28,31,30,31,30,31,31,30,31,30,31]

u_vel_stor = np.ma.zeros((122,50,70))
z_val_stor = np.ma.zeros((122,50,70))

for yr in np.arange(1996,1999+1):
    dir = '/t1/scratch/liz/tmpdir_VIP-LD.HCo13T/outputs/'+str(yr)+'/' 
    n=0
    for mon in np.arange(5,7+1):
        monstr = str(mon+1).zfill(2)
        for day in np.arange(1,mon_days[mon]+1):
            daystr = str(day).zfill(2)
            file = dir + 'VIP-LD.HCo13T_avg_' + str(yr) + '-' + monstr + '-' + daystr + 'T00:00:00.nc'
            u_vel_stor[n,:,:],z_val_stor[n,:,:] = get_vels(file)
            n+=1
    avg_u = np.mean(u_vel_stor,axis=0)
    avg_zval = np.mean(z_val_stor,axis=0)
    X = np.arange(avg_u.shape[1]) 
    plt.figure()
    plt.pcolor(X,avg_zval,avg_u,vmin=-1,vmax=1)
    plt.colorbar()
    plt.title(str(yr))
    
plt.show()
