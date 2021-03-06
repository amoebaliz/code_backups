import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
from numpy import genfromtxt
import pyroms

def convert_time(times):
     date = np.zeros(len(times))
     for nt in range(len(times)):
        newdate = ref + dt.timedelta(seconds=times[nt])
        date[nt] = pltd.date2num(newdate)
     return date
    
def get_station_temp(j):
#    sta_fil = sta_dir + 'VIP_station_' + str(stations[j]).zfill(2) + '.nc'
    sta_fil = '/t1/scratch/liz/tmpdir_VIP-LD.HCo11T/VIP-LD.HCo11T_sta.nc'
    fid = nc.Dataset(sta_fil)
    time = fid.variables['ocean_time'][:]
    temp = np.transpose(np.squeeze(fid.variables['temp'][:,19,:]))
#    u = np.transpose(np.squeeze(fid.variables['u'][:,19,:]))
    return time, temp

def plot_mld(temp, sta_depth,ax):
    dep_store = np.zeros(temp.shape[1])
    mld_temp_store = np.zeros(temp.shape[1])
    for nt in range(temp.shape[1]):
        sst = np.squeeze(temp[-1,nt])
        temp_n = np.squeeze(temp[::-1,nt])
        mld_temp = sst - 0.2
        depth = sta_depth[::-1]

        dep_store[nt] = depth[np.argmax(temp_n<mld_temp)]
        mld_temp_store[nt] = np.mean(temp_n[np.where(temp_n >= mld_temp)])

    ax.plot(plot_time,dep_store,'0.2')

    plt.figure()
    plt.plot(plot_time,mld_temp_store)

def plot_stations():

    plt.figure()
    plt.pcolor(bathy_mask)
    plt.colorbar()
    plt.plot(x_vals[stations],y_vals[stations],'o',mfc='y',mec='k',mew=1.5)
    plt.xlim(0,bathy_mask.shape[1])
    plt.ylim(0,bathy_mask.shape[0])
    plt.show()

###########
vip_grd = pyroms.grid.get_ROMS_grid('VIP')
depth = vip_grd.vgrid.z_r[0]
bathy = vip_grd.vgrid.h
mask_rho = vip_grd.hgrid.mask_rho
bathy_mask = np.ma.masked_where(mask_rho == 0, bathy)

sta_dir = '/t3/workdir/liz/MODELS/VIP/PostProc/tmpdir_VIP-LD.HCo10T/stations/'
sta_fil = '/home/liz/ANALYSES/VIP/VIP_stations.txt'
my_data=genfromtxt(sta_fil, usecols=np.arange(0,2),skip_header=2)

max_sig = 20

x_vals = my_data[:,0]
y_vals = my_data[:,1]

#stations = [4,5,6,7,8,16,19]
stations = [19]

ref = dt.datetime(1900,1,1,0,0)
levels = np.linspace(12,32,num=50)
#levels = np.linspace(-2,2,num=50)
fig,ax = plt.subplots(len(stations),sharex='col')
fig.tight_layout(h_pad = -.5)
plt.subplots_adjust(right=0.90)
cbar_ax = fig.add_axes([0.93,0.12,0.02,.75])

for nt in range(len(stations)):
    time, temp = get_station_temp(nt)
    plot_time = convert_time(time)
    sta_depth = depth[:,y_vals[stations[nt]],x_vals[stations[nt]]]

#    C = ax[nt].contourf(plot_time,sta_depth[max_sig:],temp[max_sig:],levels, hold='on')
#    ax[nt].contour(plot_time, sta_depth[max_sig:], temp[max_sig:,:],(27,),colors = 'k', linewidths = 1)
#    ax[nt].xaxis_date()
#    ax[nt].set_ylim(-100,0)
#    ax[nt].set_yticks([-100,-50,0]) 

    C = ax.contourf(plot_time,sta_depth[max_sig:],temp[max_sig:],levels, hold='on',alpha=0.7)
#    C = ax.contourf(plot_time,sta_depth[max_sig:],u[max_sig:],levels,cmap='bwr')
#    ax.contour(plot_time, sta_depth[max_sig:], temp[max_sig:,:],(27,),colors = 'k', linewidths = 1)
#    ax.contour(plot_time, sta_depth[max_sig:], u[max_sig:,:],(0,),colors = 'k', linewidth=1)
    plot_mld(temp,sta_depth,ax)
    ax.xaxis_date()
    ax.set_ylim(-100,0)
    ax.set_yticks([-100,-50,0]) 


cbar = plt.colorbar(C,cax=cbar_ax)

#plot_stations()

plt.show()



