import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.animation as animation
import datetime as dt
from numpy import genfromtxt

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(seconds=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

# SET UP
tport_file = 'vip_u_ALLm_transport_col280_rows145-178.nc'
time_vals = ['ocean_time']
color_vals = ['-b','-g']

place_time = convert_time(np.arange(35071.25,36527,7))
# INITIAL FIGURE
fig, ax = plt.subplots(1, figsize=(19,12))

# SST PLOT

# TRANSPORT PLOT
fidT = nc.Dataset(tport_file)
time = fidT.variables['ocean_time'][:]
plot_time = convert_time(time)
pos_trans = fidT.variables['pos_transport'][:]
neg_trans = fidT.variables['neg_transport'][:]
fSv = 1000000

ax.plot(plot_time,pos_trans/fSv)
ax.plot(plot_time,neg_trans/fSv)

ax.xaxis_date()
ax.xaxis.set_ticks_position('bottom')
ax.plot(((pltd.date2num(dt.datetime(1996,1,1)),pltd.date2num(dt.datetime(1999,12,31)))),(0,0),'-k')
ax.set_xlim(pltd.date2num(dt.datetime(1996,1,1)),pltd.date2num(dt.datetime(1999,12,31)))
plt.show()
