import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
from numpy import genfromtxt


my_data=genfromtxt('data.txt', usecols=np.arange(0,2))
my_data = my_data[542:589+1,:]

n=0
soi_date = np.zeros(4*12)
for nyr in range(1996,1999+1):
    for nmon in range(1,12+1):
        date_val = dt.datetime(nyr,nmon,15)
        soi_date[n] = pltd.date2num(date_val)
        n+=1

ref = dt.datetime(1900,1,1,0,0)
fid = nc.Dataset('SST_1996-1999.nc','r')
time = fid.variables['ocean_time'][:]

sst = fid.variables['temp'][:]
sst2 = np.squeeze(np.mean(np.mean(sst,axis=3),axis=2))

date2 = np.zeros(len(time))
for nt in range(len(time)):
    day_time = ref + dt.timedelta(seconds=time[nt])
    date2[nt] = pltd.date2num(day_time)

fig,ax1 = plt.subplots(1, figsize = (16,5))
ax1.plot(date2,sst2,'-b')

ax1.set_xlim(date2[0],date2[-1])
ax1.xaxis_date()
ax1.set_yticks([26,28,30])
ax1.set_ylim(18,32)
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax1.set_ylabel('oC',color='b')
ax1.yaxis.set_label_coords(-.05,.75)

soi = my_data[:,1]
la_soi = np.ma.masked_where(soi<0,soi)
el_soi = np.ma.masked_where(soi>0,soi)

ax2 = ax1.twinx()

#bar_width=(soi_date[2]-soi_date[1])
bar_width=30

ax2.bar(soi_date,el_soi,bar_width,color='r')
ax2.bar(soi_date,la_soi,bar_width,color='b')
#ax2.plot(soi_date,my_data[:,1],'-k')
ax2.plot([date2[0],date2[-1]],[0,0],'--k')
ax2.set_yticks([-2,0,2])
ax2.set_ylim(-3,8)
ax2.set_ylabel('SOI',rotation=270)
ax2.yaxis.set_label_coords(1.05,.3)

plt.title('SOI and Average VIP SST')
plt.show()

