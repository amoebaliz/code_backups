import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as pltd 

def get_time(j):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    ncfile = srcname[j] + '_VIPbound_sst_ts_A.nc'
    fid = nc.Dataset(ncfile)
    time = fid.variables[time_vals[j]][:]
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(days=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

def get_sst(j,k):
    ncfile = srcname[j] + '_VIPbound_sst_ts_' + boundID[k] + '.nc'
    fid = nc.Dataset(ncfile)
    sst = np.squeeze(fid.variables[sst_vals[j]][:])
    sst = np.ma.masked_where(sst>500,sst)
    avg_sst = np.squeeze(np.mean(sst,axis=1))
    return avg_sst

srcname = ['ct','cortad']
boundID = ['A','B','C1','C2','D']
sst_vals = ['temp','sst']
time_vals = ['ocean_time','time']

# CREATE FIGURE
fig,ax = plt.subplots(5,1,sharex = 'col')
plt.tight_layout()
# ASSEMBLE FIGURE
for nsrc in range(len(srcname)):
    plot_time = get_time(nsrc)
    for nbound in range(len(boundID)):
        sst = get_sst(nsrc,nbound)
        ax[nbound].plot(plot_time,sst)
        ax[nbound].xaxis_date()
        ax[nbound].set_ylim(22,32)
        ax[nbound].yaxis.set_ticks([22,27,32])
        ax[nbound].xaxis.set_ticks_position('bottom')
plt.show()
