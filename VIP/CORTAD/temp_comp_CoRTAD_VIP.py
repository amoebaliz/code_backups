#
import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as pltd
import datetime as dt

months = pltd.MonthLocator()
monthsFmt = pltd.DateFormatter('%m')

hrly_file = '/t3/workdir/liz/VIP/Runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150323/hrly_swrad_sst.nc'
daily_file = '/t3/workdir/liz/VIP/Runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150401/daily_swrad_sst.nc'
CoRTAD_file = '/t3/workdir/liz/scripts/CORTAD_to_VIP_remap/remapped_CORTADv4_FilledSST_to_VIP.nc'
#CTroms_file = '/t3/workdir/liz/VIP/VIP_grid/VIP_interpolations/ctroms_2007_sst.nc'


fid_3hr = nc.Dataset(hrly_file,'r')
fid_cor = nc.Dataset(CoRTAD_file,'r')
#fid_ctrom = nc.Dataset(CTroms_file,'r')
fid_daily = nc.Dataset(daily_file,'r')

ref = dt.datetime(1900,1,1,0,0)
## GET THE TIME VARIABLES FOR EACH DATA SET
time_hrly = fid_3hr.variables['ocean_time'][:]
time_daily = fid_daily.variables['ocean_time'][:]
time_cor = fid_cor.variables['time'][:]

tm_hrly = time_hrly[:]
tm_daily = time_daily[:]

tidx = np.where((time_cor > tm_hrly[0]-1) & (time_cor < tm_hrly[-1]+1))
tm_cor = time_cor[tidx[0]]
## GET THE TEMP VALUES FOR EACH DATASET
sst_hrly = np.squeeze(fid_3hr.variables['temp'][:])
sst_daily = np.squeeze(fid_daily.variables['temp'][:])
sst_cor = fid_cor.variables['sst'][tidx[0],:,:]
sst_cor[sst_cor<-10] = float('NaN')

#sst_vip_mean = np.mean(np.mean(sst_vip,axis=2),axis=1)
#sst_ctrom_mean = np.mean(np.mean(sst_ctroms,axis=2),axis=1)
#sst_cor_mean = np.nanmean(np.nanmean(sst_cor,axis=2),axis=1)
#temp_dif = sst_vip_mean - sst_cor_mean
#temp_dif_std = np.std(temp_dif)
#temp_dif_mean = np.mean(temp_dif) 
#print temp_dif_mean
#print temp_dif_std

dates1 = [pltd.date2num(ref + dt.timedelta(days=tm_hrly[nt])) for nt in range(len(tm_hrly))]
dates2 = [pltd.date2num(ref + dt.timedelta(days=tm_daily[nt])) for nt in range(len(tm_daily))]
dates3 = [pltd.date2num(ref + dt.timedelta(days=tm_cor[nt])) for nt in range(len(tm_cor))]

#fig,ax = plt.subplots()

#ax.plot_date(dates1,sst_vip_mean,'-k')
#ax.plot_date(dates2,sst_ctrom_mean,'-r')
#ax.plot_date(dates3,sst_cor_mean,'-g')
#plt.title('Weekly SST Averaged Over VIP Domain')
#plt.ylabel('Average Weekly SST ($^o$C)')
#plt.ylim(24,38)
#plt.legend(['VIP','CT-ROMS','CoRTAD'])

#ax.xaxis.set_major_locator(months)
#ax.xaxis.set_major_formatter(monthsFmt)
#fig.autofmt_xdate()
#ax.fmt_xdata = pltd.DateFormatter('%Month')

##
## SUBSET OF VIP DOMATIN COMPARISON ##
##

#x1 = 50; x2 = 80 
#y1 = 140; y2 = 160

#x1 = 550; x2 = 570 
#y1 = 90; y2 = 110

#ix = np.arange(x1,x2+1)
#iy = np.arange(y1,y2+1)

#fid_grid = nc.Dataset('/t3/workdir/liz/VIP/VIP_grid/VIP_grd_high_res_bathy.nc')
#mask = fid_grid.variables['mask_rho'][:]

#sst_vip_mean = np.mean(np.mean(sst_vip[:,y1:y2,x1:x2],axis=2),axis=1)
#sst_ctrom_mean = np.mean(np.mean(sst_ctroms[:,y1:y2,x1:x2],axis=2),axis=1)
#sst_cor_mean = np.nanmean(np.nanmean(sst_cor[:,y1:y2,x1:x2],axis=2),axis=1)

sst_hrly_mean = np.mean(np.mean(sst_hrly,axis=2),axis=1)
sst_daily_mean = np.mean(np.mean(sst_daily,axis=2),axis=1)
sst_cor_mean = np.nanmean(np.nanmean(sst_cor,axis=2),axis=1)

fig, ax = plt.subplots(1, figsize=(16,4))
#plt.subplot(121)
#plt.pcolor(mask)

#plt.plot([x1,x2],[y1,y1],'w')
#plt.plot([x1,x2],[y2,y2],'w')
#plt.plot([x1,x1],[y1,y2],'w')
#plt.plot([x2,x2],[y1,y2],'w')
#plt.xlim(0,602)
#plt.ylim(0,342)

#plt.subplot(122)
ax.plot_date(dates1,sst_hrly_mean,'-r',label='VIP_3hrly_swrad')
ax.plot_date(dates2,sst_daily_mean,'-k',label='VIP_daily_swrad')
ax.plot_date(dates3,sst_cor_mean,'-g',linewidth=2.0,label='CoRTAD')
ax.yaxis.tick_left()
ax.xaxis.tick_bottom()
ax.yaxis.set_ticks([24,26,28,30,32])

plt.ylim(24,32)
plt.xlim(np.min(dates1),np.max(dates1))
plt.title('SST Averaged Over VIP Domain')
plt.ylabel('Average SST ($^o$C)')
plt.legend(loc=2)
fig.autofmt_xdate()
#ax.xaxis.set_major_locator(months)
#ax.xaxis.set_major_formatter(monthsFmt)
#fig.autofmt_xdate()



plt.show()
