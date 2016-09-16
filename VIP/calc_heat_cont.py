# Calculate the heat content above a particular depth level
import pyroms as py
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as pltd

fid = nc.Dataset('/t3/workdir/liz/VIP/VIP_grid/VIP_grd_high_res_bathy.nc','r')
h = fid.variables['h'][:]; h=h[:]
pn = fid.variables['pn'][:]; dx = 1/pn[:]
pm = fid.variables['pm'][:]; dy = 1/pm[:]
fid.close()
ref = dt.datetime(1900,1,1,0,0)

n_dep = 100
n_dep = n_dep*-1
ct_num = np.arange(17900,18262+1)
#ct_num = np.arange(17900,17942+1)
heat_cont_ctr_up = np.zeros(len(ct_num)) 
heat_cont_ctr_low = np.zeros(len(ct_num))
time_ctr = np.zeros(len(ct_num)) 
n = 0

ctrom_dir = '/t3/workdir/liz/VIP/VIP_grid/VIP_interpolations/ctroms_vel/' 

for num in ct_num:
   ncfile = ctrom_dir + 'coral_avg_' + str(num) + '.nc'
   fid = nc.Dataset(ncfile,'r')
   hc = fid.variables['hc'][:]; hc=hc[:]
   s_w = fid.variables['s_w'][:]; s_w=s_w[:]
   Cs_w = fid.variables['Cs_w'][:]; Cs_w=Cs_w[:]
   zeta = fid.variables['zeta'][:]; zeta=zeta[:]
   time = fid.variables['ocean_time'][:]; tmp=time[:]
   temp = fid.variables['temp'][:]; temp = np.squeeze(temp[:])
   
   z_w = py.vgrid.z_w(h, hc, 51, s_w, Cs_w, zeta, 4); z_w = z_w[:] 
   # z_w[z_w>100]=np.nan 
   # NOTE: nans for land not necessary b/c temp values are masked AND 
   # all z_w values are the same so dz = 0.

   # for ABOVE a certain value use greater than (>)
   # READS: all values greater than n_dep, keep original z_w values. All others set to n_dep
   # reverse this convention for heat content below a certain depth
   z_w_1 = np.where(z_w>n_dep,z_w,n_dep)
   z_w_2 = np.where(z_w<n_dep,z_w,n_dep)

   dz_up = z_w_1[1:,:,:]-z_w_1[:-1,:,:]   
   dz_low = z_w_2[1:,:,:]-z_w_2[:-1,:,:]

   vol_up = dz_up*dx*dy
   vol_low = dz_low*dx*dy

   heat_cont_up = np.sum(np.sum(np.sum(vol_up*temp,axis=2),axis=1),axis=0)
   total_vol_up = np.sum(np.sum(np.sum(vol_up,axis =2),axis=1),axis=0)
   heat_cont_low = np.sum(np.sum(np.sum(vol_low*temp,axis=2),axis=1),axis=0)
   total_vol_low = np.sum(np.sum(np.sum(vol_low,axis =2),axis=1),axis=0)   

   newdate = pltd.date2num(ref + dt.timedelta(seconds=tmp[0]))
   time_ctr[n] = newdate 

   heat_cont_ctr_up[n] = heat_cont_up/total_vol_up
   heat_cont_ctr_low[n] = heat_cont_low/total_vol_low
   n = n+1
   fid.close()


#vip_dir = '/t3/workdir/liz/VIP/VIP_runs/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150205/'
# vip_num = np.arange(5,7+1)
#vip_num = np.arange(5,78+1)
#heat_cont_vip_up = np.zeros(len(vip_num))
#heat_cont_vip_low = np.zeros(len(vip_num))
#time_vip = np.zeros(len(vip_num))
#m = 0

#for num in vip_num:

#   ncfile = vip_dir + 'VIP-LD.HCo07T_avg_' + str(num).zfill(4) + '.nc'
#   fid = nc.Dataset(ncfile,'r')
#   hc = fid.variables['hc'][:]; hc=hc[:]
#   s_w = fid.variables['s_w'][:]; s_w=s_w[:]
#   Cs_w = fid.variables['Cs_w'][:]; Cs_w=Cs_w[:]
#   zeta = fid.variables['zeta'][:]; zeta=zeta[:]
#   time = fid.variables['ocean_time'][:]; tmp=time[:]
#   temp = fid.variables['temp'][:]; temp = np.squeeze(temp[:])

#   z_w = py.vgrid.z_w(h, hc, 51, s_w, Cs_w, zeta, 4); z_w = z_w[:] 
#   z_w_1 = np.where(z_w>n_dep,z_w,n_dep)
#   z_w_2 = np.where(z_w<n_dep,z_w,n_dep)
   
#   dz_up = z_w_1[1:,:,:]-z_w_1[:-1,:,:]
#   dz_low = z_w_2[1:,:,:]-z_w_2[:-1,:,:]
   
#   vol_up = dz_up*dx*dy
#   vol_low = dz_low*dx*dy

#   heat_cont_up = np.sum(np.sum(np.sum(vol_up*temp,axis=2),axis=1),axis=0)
#   total_vol_up = np.sum(np.sum(np.sum(vol_up,axis =2),axis=1),axis=0)
#   heat_cont_low = np.sum(np.sum(np.sum(vol_low*temp,axis=2),axis=1),axis=0)
#   total_vol_low = np.sum(np.sum(np.sum(vol_low,axis =2),axis=1),axis=0)

#   newdate = pltd.date2num(ref + dt.timedelta(seconds=tmp[0]))
#   time_vip[m] = newdate

#   heat_cont_vip_up[m] = heat_cont_up/total_vol_up
#   heat_cont_vip_low[m] = heat_cont_low/total_vol_low
   
#   m = m+1
#   fid.close()

vip_dir = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150313/'

#vip_num = np.arange(5,5+1)
vip_num = np.arange(5,304+1)
heat_cont_vip_up = np.zeros(len(vip_num))
heat_cont_vip_low = np.zeros(len(vip_num))
time_vip = np.zeros(len(vip_num))
o = 0

for num in vip_num:

   ncfile = vip_dir + 'VIP-LD.HCo07T_avg_' + str(num).zfill(4) + '.nc'
   fid = nc.Dataset(ncfile,'r')
   hc = fid.variables['hc'][:]; hc=hc[:]
   s_w = fid.variables['s_w'][:]; s_w=s_w[:]
   Cs_w = fid.variables['Cs_w'][:]; Cs_w=Cs_w[:]
   zeta = fid.variables['zeta'][:]; zeta=zeta[:]
   time = fid.variables['ocean_time'][:]; tmp=time[:]
   temp = fid.variables['temp'][:]; temp = np.squeeze(temp[:])

   z_w = py.vgrid.z_w(h, hc, 51, s_w, Cs_w, zeta, 4); z_w = z_w[:]
   z_w_1 = np.where(z_w>n_dep,z_w,n_dep)
   z_w_2 = np.where(z_w<n_dep,z_w,n_dep)
   
   dz_up = z_w_1[1:,:,:]-z_w_1[:-1,:,:]
   dz_low = z_w_2[1:,:,:]-z_w_2[:-1,:,:]
   
   vol_up = dz_up*dx*dy
   vol_low = dz_low*dx*dy

   heat_cont_up = np.sum(np.sum(np.sum(vol_up*temp,axis=2),axis=1),axis=0)
   total_vol_up = np.sum(np.sum(np.sum(vol_up,axis =2),axis=1),axis=0)
   heat_cont_low = np.sum(np.sum(np.sum(vol_low*temp,axis=2),axis=1),axis=0)
   total_vol_low = np.sum(np.sum(np.sum(vol_low,axis =2),axis=1),axis=0)

   newdate = pltd.date2num(ref + dt.timedelta(seconds=tmp[0]))
   time_vip[o] = newdate

   heat_cont_vip_up[o] = heat_cont_up/total_vol_up
   heat_cont_vip_low[o] = heat_cont_low/total_vol_low

   o = o+1
   fid.close()

fig, axarr = plt.subplots(2,sharex=True)
a = axarr[0].plot_date(time_ctr,heat_cont_ctr_up,'-k',label='CT-ROMS')
b = axarr[0].plot_date(time_vip,heat_cont_vip_up,'-r',label='VIP-old')
axarr[0].set_ylabel('Temperature($^o$C)')
title = 'Average Temperature a) Above and b) Below ' + str(-1*n_dep) + ' Meters'
axarr[0].set_title(title)

axarr[1].plot_date(time_ctr,heat_cont_ctr_low,'-k',label='CT-ROMS')
axarr[1].plot_date(time_vip,heat_cont_vip_low,'-r',label='VIP-old')
axarr[1].set_ylabel('Temperature($^o$C)')
#axarr[1].legend(loc=3)
y = np.linspace(10,11,5,endpoint=True)
fig.autofmt_xdate()

fig.savefig('heat_content_month1.eps',format='eps', dpi=1000)
plt.show()

