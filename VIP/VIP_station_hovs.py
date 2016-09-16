import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import pyroms as py
from numpy import genfromtxt

nsta = 19

sta_dir = '/home/liz/ANALYSES/VIP/'
sta_fil = sta_dir + 'VIP_stations.txt'

my_data=genfromtxt(sta_fil, usecols=np.arange(0,2),skip_header=2)
fid = nc.Dataset('/t1/scratch/liz/Inputs/VIP-LD.HCo11T/Grid/VIP_grd_high_res_bathy_interp2.nc','r')
#h = fid.variables['h'][:]
#mask = fid.variables['mask_rho'][:]
#bathy = np.ma.zeros(h.shape[0].h.shape[1])
#bathy = np.ma.masked_where(mask == 0, h)

#x_vals = my_data[:,0]
#y_vals = my_data[:,1]

fid_sta = nc.Dataset('/t1/scratch/liz/tmpdir_VIP-LD.HCo11T/VIP-LD.HCo11T_sta.nc','r')
temp = np.squeeze(fid_sta.variables['temp'][:,nsta,:])
time_val = fid_sta.variables['ocean_time'][:]

#Cs_r = fid_sta.variables['Cs_r'][:]
#zeta = fid_sta.variables['zeta'][:]
#VIP = py.grid.get_ROMS_grid('VIP')

#n_sta = len(x_vals)
ref = dt.datetime(1900,1,1,0,0)

def plot_date(times):
    date = np.zeros(len(times))
    for nt in range(len(times)):
        newdate = ref + dt.timedelta(seconds=times[nt])
        date[nt] = pltd.date2num(newdate)
    return date

plt_dates = plot_date(time_val)
s_rho = np.arange(0,50)
x = plt_dates
X, Y = np.meshgrid(x,s_rho)

for n in np.arange(19):
    temp_sta = np.squeeze(temp[:])
    fig, ax = plt.subplots(1)
    p = ax.pcolor(X,Y,temp_sta.T)
    plt.colorbar(p)
    plt.title('Temperature oC')

    ax.xaxis.set_major_formatter(pltd.DateFormatter('%m/%d/%y'))
    ax.xaxis.set_major_locator(pltd.DayLocator())
    plt.gcf().autofmt_xdate()

    ax.set_xlim(plt_dates[0],plt_dates[-1])
    ax.set_xticks(np.arange(plt_dates[0],plt_dates[-1],50))
    ax.set_ylim(0,s_rho[-1])
    ax.set_ylabel('s_rho')

    # Map station 
#    plt.figure()
#    plt.pcolor(bathy)
#    plt.plot(x_vals[n],y_vals[n],'o',mfc='y',mec='k',mew=1.5)
#    plt.xlim(0,h.shape[1])
#    plt.ylim(0,h.shape[0])
#    plt.colorbar()
fig.savefig('hov_test.eps',format='eps', dpi=1000)
#    plt.title('Station #' +str(n+1))
