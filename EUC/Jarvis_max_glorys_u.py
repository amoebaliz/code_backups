import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
import csv

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN SECONDS SINCES (1991,12,4,0,0)
    ref = dt.datetime(1991,12,04,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(seconds=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals


fid = nc.Dataset('/t3/workdir/liz/external_data/GLORYS_EUC/GLORYS_jarvis_92_13_U.nc','r')

time = fid.variables['time_counter'][:]
depth = depth = fid.variables['deptht'][:]
Idep = np.where((10<depth)&(depth<300))

u = fid.variables['vozocrtx'][:,Idep[0],:,:]
time = fid.variables['time_counter'][:]

plot_dates = convert_time(time)

u_avg = np.squeeze(np.mean(u,axis=2))
u_max = np.max(u_avg,axis=1)

#u_max_ma = np.ma.masked_where(u_max>10,u_max)

plot_dates = convert_time(time)

with open('jarvis_max_u.csv', 'wb') as f:
     a = csv.writer(f)
     for nt in range(len(u_max)):
         a.writerow([float(u_max[nt])])

fig, ax = plt.subplots()
ax.plot(plot_dates,u_max)
ax.xaxis_date()
plt.show()

