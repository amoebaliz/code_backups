import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(days=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals


file = '/data/external/P1/Data/CORTAD/Version4/cortadv4_TSA_DHW_coral_1998-1999.nc'
fid = nc.Dataset(file)

time = fid.variables['time'][:60]
plot_time = convert_time(time)

DHW = fid.variables['dhw'][:60,:,:]
DHW = np.ma.masked_where(DHW<0,DHW)

DHW_max = np.max(DHW,axis=0)
DHW_max_when = np.argmax(DHW,axis=0)

when_dates = plot_time[DHW_max_when]

plt.figure()
plt.pcolor(DHW_max,vmin=0,vmax=16)
plt.colorbar()

plt.figure()
plt.pcolor(when_dates)
plt.colorbar()
plt.show()

