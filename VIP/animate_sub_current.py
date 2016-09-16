import netCDF4 as nc
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime as dt
import pyroms

def get_vels(ncfile):
     fid = nc.Dataset(ncfile,'r')
     u = np.squeeze(fid.variables['u'][:,:,130:200,280])
     zeta = np.squeeze(fid.variables['zeta'][:]) 
     z = pyroms.vgrid.z_w(VIP.vgrid.h, VIP.vgrid.hc, VIP.vgrid.N, VIP.vgrid.s_rho, VIP.vgrid.Cs_r, zeta, VIP.vgrid.Vtrans)
     z_vals = np.squeeze(z[:,:,130:200,290])
     z_vals = np.ma.masked_where(z_vals>0,z_vals)
     #return u_vel_2, v_vel_2, cur_speed
     return u, z_vals

VIP = pyroms.grid.get_ROMS_grid('VIP')

mon_days = [31,28,31,30,31,30,31,31,30,31,30,31]

for yr in np.arange(1997,1997+1):
    dir = '/data/external/P6/ROMS/VIP/VIP-LD.HCo12T/1997/'
    #dir = '/t1/scratch/liz/tmpdir_VIP-LD.HCo13T/outputs/'+str(yr)+'/'
    file0 = 'VIP-LD.HCo12T_his_' + str(yr) + '-04-01T00:00:00.nc'
    file = dir+file0
    print file
    u_vel,z_vals = get_vels(file)
    X = np.arange(u_vel.shape[1])

    # FIGURE SETUP
    fig = plt.figure(figsize=(20,10))
    ax1 = plt.subplot(axisbg=[.2,.2,.2])
    plt.tight_layout(pad=2.0, w_pad=1, h_pad=2.0)
    pcol = ax1.pcolor(X,z_vals,u_vel,vmin=20,vmax=30)
    cbar = plt.colorbar(pcol)

    ims1 = []
    #if yr%4 == 0:
    #   mon_days[1] = 29
    #else:
    #   mon_days[1] = 28 
    for mon in np.arange(3,4):
#    for mon in np.arange(1):
        monstr = str(mon+1).zfill(2)
        for day in np.arange(1,mon_days[mon]+1):   
            print day
            daystr = str(day).zfill(2)
            for hr in range(23):
                hrstr = str(hr).zfill(2)
                file_2 = dir + 'VIP-LD.HCo12T_his_' + str(yr) + '-' + monstr + '-' + daystr + 'T' + hrstr + ':00:00.nc'
                print file_2
                u_vel_2,z_vals_2 = get_vels(file_2)
                ims1.append((ax1.pcolor(X,z_vals_2,u_vel_2,vmin=20,vmax=30),))
    gif_file = 'VIP_' + str(yr) + '_280_sub_hrly_sub.gif'
    im_ani = animation.ArtistAnimation(fig, ims1, blit=True)
    im_ani.save(gif_file, writer = 'imagemagick',fps=10)
#plt.show()
