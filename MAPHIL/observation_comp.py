import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.animation as animation
import datetime as dt
from numpy import genfromtxt

def convert_time(time,ref):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(seconds=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

logg_data=genfromtxt('/t3/workdir/liz/external_data/MaPhil_Temp_Loggers/Visca1_2013_05_22.txt',delimiter=",", skip_header=2)
logg_sst = logg_data[:,2]
logg_incr = 30*60*(np.arange(0,logg_data.shape[0]))
ref = dt.datetime(2012,5,7,0,0)
plot_logg_dates = convert_time(logg_incr,ref)
end_val = -3

logg_data2=genfromtxt('/t3/workdir/liz/external_data/MaPhil_Temp_Loggers/Visca3_2015_01_06.txt',delimiter=",", skip_header=2)
logg_sst2 = logg_data[:,2]
logg_incr2 = 30*60*(np.arange(0,logg_data.shape[0]))
ref = dt.datetime(2014,6,26,12,0)
plot_logg_dates2 = convert_time(logg_incr,ref)

sta_file = '/t1/scratch/liz/tmpdir_MaPhil-LD.HCo05T/MaPhil-LD.HCo05T_sta.nc'

#avg_file = 'MaPhil_sst2.nc'
fid = nc.Dataset(sta_file)
sta_sst = np.squeeze(fid.variables['temp'][5:,-1,-1])
sta_time = fid.variables['ocean_time'][5:]
plot_sta_dates = convert_time(sta_time,dt.datetime(1900,1,1,0,0))

# FIGURE
fig,ax = plt.subplots(1)
ax.plot(plot_logg_dates[:end_val],logg_sst[:end_val])
ax.plot(plot_logg_dates2[:end_val],logg_sst2[:end_val],'-b')
ax.plot(plot_sta_dates,sta_sst,'-r')
ax.xaxis_date()
plt.show()    

