import netCDF4 as nc
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime as dt

def get_vels(file):
     fid = nc.Dataset(file,'r')
     u_vel = np.squeeze(fid.variables['u'][:,-1,:,:])
     v_vel = np.squeeze(fid.variables['v'][:,-1,:,:])
     # Interpolate
     u_vel_2 = (u_vel[1:-1,1:] + u_vel[1:-1,:-1])/2
     v_vel_2 = (v_vel[1:,1:-1] + v_vel[:-1,1:-1])/2
     # Speed Calculation
     cur_speed = np.sqrt((u_vel_2)**2+(v_vel_2)**2)
     return u_vel_2, v_vel_2, cur_speed


# FILE DETAILS
dir = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo11T/'

# Initial image file
file0 = '1995/VIP-LD.HCo11T_avg_1995-12-31T00:00:00.nc'
file = dir+file0

u_vel, v_vel,cur_speed = get_vels(file)

# FIGURE SETUP
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot(axisbg=[.2,.2,.2])
ax1.set_ylim(0,u_vel.shape[0])
ax1.set_xlim(0,u_vel.shape[1])
plt.tight_layout(pad=2.0, w_pad=1, h_pad=2.0)

# INITIAL FRAME
pcol = ax1.pcolor(cur_speed,vmin=0,vmax=1.0)
cbar = plt.colorbar(pcol)
cbar.ax.set_ylabel('Speed (m/s)', fontsize =14, labelpad = 20, rotation=270)

print np.min(cur_speed), np.max(cur_speed)

x = np.arange(0,cur_speed.shape[1])
y = np.arange(0,cur_speed.shape[0])
X, Y = np.meshgrid(x, y)
afreq = 10
Q = ax1.quiver(X[::afreq,::afreq],Y[::afreq,::afreq],u_vel[::afreq,::afreq],v_vel[::afreq,::afreq],width = 0.0015)

#plt.show()
ims1 = []

mon_days = [31,28,31,30,31,30,31,31,30,31,30,31]

for yr in np.arange(1996,1996+1):
    yrstr = str(yr)

    for mon in np.arange(9,12+1):
#    for mon in np.arange(1,2):
        monstr = str(mon).zfill(2)

        for day in np.arange(1,mon_days[mon-1]+1):   
#        for day in np.arange(1,2):   
            daystr = str(day).zfill(2)
 
            file = dir + yrstr + '/' + 'VIP-LD.HCo11T_avg_' + yrstr + '-' + monstr + '-' + daystr + 'T00:00:00.nc' 
            print file
            u_vel, v_vel, cur_speed = get_vels(file)
 
            ims1.append((ax1.pcolor(cur_speed,vmin=0,vmax=1),\
                         ax1.quiver(X[::afreq,::afreq],Y[::afreq,::afreq],u_vel[::afreq,::afreq],v_vel[::afreq,::afreq],width = 0.0015),\
                         ))

im_ani = animation.ArtistAnimation(fig, ims1, blit=True)
im_ani.save('VIP_currents.gif', writer = 'imagemagick',fps=10)
#plt.show()
