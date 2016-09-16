import numpy as np
import netCDF4 as nc
from numpy import genfromtxt
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd

my_data=genfromtxt('/t3/workdir/liz/external_data/NOAA_SOI/soi.txt', usecols=np.arange(1,13))

# yrs = np.arange(1951,2016+1)
# yr = 1983
# yr_I = np.where(yrs == yr)

bar_widths = [31,28,31,30,31,30,31,31,30,31,30,31]
strt_date = pltd.date2num(dt.datetime(1951,1,1))
end_date = pltd.date2num(dt.datetime(2015,12,31))

# my_data = my_data[yr_I[0],:]
soi = my_data.flatten()
la_soi = np.ma.masked_where(soi<0,soi)
el_soi = np.ma.masked_where(soi>0,soi)

fig, ax2 = plt.subplots(1, figsize=(20,3))

#soi_date = np.zeros(12)
#bar_widths2 = np.zeros(12)
#check_dates = np.zeros(12)
n=0
for nyr in range(1951,2015+1):
    if nyr%4 == 0:
       bar_widths[1] = 29
    else:
       bar_widths[1] = 28

    for nmon in range(12):
        date_val = dt.datetime(nyr,nmon+1,1) 
        soi_date = pltd.date2num(date_val)
        bar_widths2 = bar_widths[nmon]
        ax2.bar(soi_date,el_soi[n],bar_widths[nmon],color='r',alpha=0.4)
        ax2.bar(soi_date,la_soi[n],bar_widths[nmon],color='b',alpha=0.4)
        n+=1
ax2.set_yticks([-2,0,2])
ax2.set_ylim(-3.6,3.6)
ax2.yaxis.set_ticks_position('left')

ax2.xaxis_date()

#ax2.set_xticks([dt.datetime(yr,1,1), dt.datetime(yr,3,1), dt.datetime(yr,5,1), dt.datetime(yr,7,1), dt.datetime(yr,9,1), dt.datetime(yr,11,1)])
#ax2.set_xticklabels(['Jan','Mar','May','Jul','Sep','Nov'])
ax2.xaxis.set_ticks_position('bottom')
ax2.set_xlim(strt_date,end_date)

plt.show()

