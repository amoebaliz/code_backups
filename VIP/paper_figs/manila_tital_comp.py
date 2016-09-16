import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd

def convert_time(times):
     ref = dt.datetime(1900,1,1,0,0,0)
     date = np.zeros(len(times))
     for nt in range(len(times)):
        newdate = ref + dt.timedelta(seconds=times[nt])
        date[nt] = pltd.date2num(newdate)
     return date

def get_hourly_times(yr_start,yr_end):
    date = []
    ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
    for nyr in range(yr_start,yr_end+1):
        if nyr%4 == 0:
           ndays[1]=29
        else:
           ndays[1]=28

        for nmon in range(12):
            # Iterates from 0 to 11 :. need +1 in datetime call 
            for nday in range(ndays[nmon]):
                # Iterates from 0 to ndays[nmon] :. need +1 in datetime call 
                for nhr in range(24):
                # Iterates from 0 to 23 :.do NOT need +1 in datetime call
                    date.append(pltd.date2num(dt.datetime(nyr,nmon+1,nday+1,nhr)))  
    return date
sta_fil ='/data/external/P6/ROMS/VIP/VIP-LD.HCo11T/stations/VIP-LD.HCo11T_sta.nc'
fid = nc.Dataset(sta_fil)

start_i = 11568
end_i = 151820

#vip_ssh = fid.variables['zeta'][:,13]
vip_ssh = fid.variables['zeta'][start_i:end_i+1:4,13]
vip_ssh = vip_ssh - np.mean(vip_ssh)
vip_time = fid.variables['ocean_time'][start_i:end_i+1:4]
plt_vip_dates = convert_time(vip_time)


tid_dir = '/data/external/P1/Data/tide_gauge/Manila/'
tid_fils = ['h370a96.dat','h370a97.dat','h370a98.dat','h370a99.dat']

arr = []
for nfil in range(len(tid_fils)):
    filename = tid_dir+tid_fils[nfil]
    arr_tmp=np.loadtxt(filename,skiprows=1,usecols=(3,4,5,6,7,8,9,10,11,12,13,14))
    arr = np.append(arr,arr_tmp.flatten())

arr = np.ma.masked_where(arr > 9000,arr)
tid_mean = np.mean(arr)

arr_norm_m = (arr - tid_mean)/1000.
plot_dates = get_hourly_times(1996,1999)

R_val = np.ma.corrcoef(arr_norm_m,vip_ssh)
print R_val[0,-1]     

 
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(plot_dates,arr_norm_m,'-r')
ax.plot(plt_vip_dates,vip_ssh, '-b')
ax.xaxis_date()
plt.show()

