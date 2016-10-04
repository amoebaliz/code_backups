import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
from matplotlib.ticker import MultipleLocator

def convert_time(times):
     #ref = dt.datetime(1900,1,3,2,0,0)
     ref = dt.datetime(1900,1,1,0,0,0)
     date = np.zeros(len(times))
     for nt in range(len(times)):
        newdate = ref + dt.timedelta(seconds=times[nt])
#        if ((nt == 25560) or (nt == 25559)): 
#           print newdate
        date[nt] = pltd.date2num(newdate)
     return date

def get_hourly_times(start_date,ndata):
    ref = start_date
    date = []
    for nt in range(ndata):
        date.append(pltd.date2num(ref + dt.timedelta(hours=nt))) # + dt.timedelta(hours=12)))

    return date

def get_manila_tides(i):
    man_tids = []
    for nt in range(i):
        tid_fil = tid_dir + tidal_files[nt]
#        if nt == 0:
#           arr_tmp=np.loadtxt(tid_fil,skiprows=489,usecols=(3,4,5,6,7,8,9,10,11,12,13,14))
#        else:
        arr_tmp=np.loadtxt(tid_fil,skiprows=1,usecols=(3,4,5,6,7,8,9,10,11,12,13,14))
        arr = arr_tmp.flatten()
        arr = np.ma.masked_where(arr > 9000,arr)
        man_tids = np.ma.append(man_tids,arr)
    tid_mean = np.mean(man_tids)
    arr_norm_m = (man_tids - tid_mean)/1000.

    return arr_norm_m

a = 11568
b = -820
#a =0
#b = -1
sta_fil ='/t1/scratch/liz/tmpdir_VIP-LD.HCo13T/VIP-LD.HCo13T_sta.nc'
fid = nc.Dataset(sta_fil)
vip_ssh = fid.variables['zeta'][a:b:4,14]
vip_ssh = vip_ssh - np.mean(vip_ssh)
vip_time = fid.variables['ocean_time'][a:b:4]
plt_vip_dates = convert_time(vip_time)


#tidal_files = ['h370a95.dat','h370a96.dat','h370a97.dat','h370a98.dat','h370a99.dat']
tidal_files = ['h370a96.dat','h370a97.dat','h370a98.dat','h370a99.dat']
tid_dir =  '/data/external/P1/Data/tide_gauge/Manila/'


arr_norm_m = get_manila_tides(len(tidal_files))
plt_dates = get_hourly_times(dt.datetime(1996,1,1,0,0,0),len(arr_norm_m))

#start_man = 60
#R_val = np.ma.corrcoef(arr_norm_m[start_man:start_man+len(vip_ssh[36:])],vip_ssh[36:])
#print R_val[0,-1]

aa = 17544 # JAN 01 00:00:00 1998
aa = 18288 # FEB 01 00:00:00 1998
#aa = 18960 # MAR 01 00:00:00 1998
aa = 19704 # APR 01 00:00:00 1998
#aa = 20424 # MAY 01 00:00:00 1998
#aa = 21168 # JUN 01 00:00:00 1998
#aa = 21888 # JUL 01 00:00:00 1998
#aa = 22632 # AUG 01 00:00:00 1998
#aa = 23376 # SEP 01 00:00:00 1998
#aa = 24096 # OCT 01 00:00:00 1998
#aa = 24840 # NOV 01 00:00:00 1998
#aa = 25560 # DEC 01 00:00:00 1998

bb = 18287 # JAN 31 23:00:00 1998
bb = 18959 # FEB 28 23:00:00 1998
#bb = 19703 # MAR 31 23:00:00 1998
bb = 20423 # APR 30 23:00:00 1998
#bb = 21167 # MAY 31 23:00:00 1998 
#bb = 21887 # JUN 30 23:00:00 1998
#bb = 22631 # JUL 31 23:00:00 1998
#bb = 23375 # AUG 31 23:00:00 1998
#bb = 24095 # SEP 30 23:00:00 1998
#bb = 24839 # OCT 31 23:00:00 1998
#bb = 25559 # NOV 30 23:00:00 1998
#bb = 26303 # DEC 31 23:00:00 1998

R_val = np.ma.corrcoef(arr_norm_m[aa:bb+1],vip_ssh[aa:bb+1])
print R_val[0,-1]

minorLocator = MultipleLocator(.25)

fig = plt.figure(figsize = (12,6))
ax = fig.add_subplot(111)


m = ax.plot(plt_dates[aa:bb+1],arr_norm_m[aa:bb+1],'-r',label='Manila gauge')
v = ax.plot(plt_vip_dates[aa:bb+1],vip_ssh[aa:bb+1], '-b', label='VIP')
#ax.plot(plt_dates[60:],arr_norm_m[60:],'-r')
#ax.plot(plt_vip_dates[36:],vip_ssh[36:], '-b')
ax.xaxis_date()
ax.set_ylim(-.8,.8)
ax.set_yticks([-.5,0,.5])
ax.set_xlim(plt_dates[aa],plt_dates[bb])
days = pltd.DayLocator()
daysFmt = pltd.DateFormatter('April %d')
ax.xaxis.set_major_locator(days)
ax.xaxis.set_major_formatter(daysFmt)
#ax.format_xdata = pltd.DateFormatter('%m/%d')
plt.tick_params(axis='both', which='both', top='off',right='off',labelsize=14)
ax.set_xticks([plt_dates[aa+24],plt_dates[aa+8*24],plt_dates[aa+15*24],plt_dates[aa+22*24],plt_dates[aa+29*24]])
ax.set_ylabel ('Water Level (m)',labelpad=20,fontsize=18)
ax.set_xlabel ('Date',labelpad=30,fontsize=18)
ax.yaxis.set_minor_locator(minorLocator)
plt.tight_layout(h_pad=2)
plt.legend(ncol=2,loc=2)
plt.savefig('tidal_comp.png')
plt.show()
