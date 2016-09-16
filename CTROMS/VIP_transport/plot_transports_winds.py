# PLOT THE WINDS AND TRANSPORTS ON A COMMON TIME AXIS

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt

wind_file = '/t3/workdir/liz/external_data/MERRA/MERRA_VIP_wind.nc'
transport_file = '/t3/workdir/liz/MODELS/CTROMS/PostProc/TRANSPORT_VIP/ctroms_u_transport_col395.nc'

mer_fid = nc.Dataset(wind_file,'r')
ct_fid = nc.Dataset(transport_file,'r')

ref = dt.datetime(1900,1,1,0,0)

def convert_time(time,ref,incr):
   date = np.zeros(len(time))
   for nt in range(len(time)):
       if incr == 'seconds':
          day_time = ref + dt.timedelta(seconds=time[nt])
       elif incr == 'days': 
          day_time = ref + dt.timedelta(days=time[nt])
       date[nt] = pltd.date2num(day_time)
   return date

Uwind = mer_fid.variables['Uwind'][:]
Vwind = mer_fid.variables['Vwind'][:]
mer_time = mer_fid.variables['time'][:]

Uwind = np.squeeze(np.mean(np.mean(Uwind,axis=2),axis=1))
Vwind = np.squeeze(np.mean(np.mean(Vwind,axis=2),axis=1))

ct_time = ct_fid.variables['ocean_time'][:]
pos_transport = ct_fid.variables['pos_transport'][:]
neg_transport = ct_fid.variables['neg_transport'][:]

ct_plt_time = convert_time(ct_time,ref,'seconds')
mer_plt_time = convert_time(mer_time,ref,'days')



def clim_vals(data,init_yr,last_yr,time):
   
   init_ref = pltd.date2num(dt.datetime(init_yr,1,1,0,0))
   last_ref = pltd.date2num(dt.datetime(last_yr+1,1,1,0,0))
 
   itime = np.where((time>=init_ref) & (time<=last_ref))
   data = data[itime[0]]
   yr = range(init_yr,last_yr+1)   

   nyr = len(yr)
   data_store = np.empty((nyr,365)); data_store[:]=np.NAN
   leap = ((init_yr%4+3)/4)
   init = 0; last_val = 365 - leap

   for nt in range(nyr):

       vals = data[init:last_val+1] 
       if leap == 0:
          vals = np.interp(np.arange(0.5,365.5),np.arange(0,366),vals)  

       # store values
       data_store[nt,0:len(vals)] = vals
       # update

       init = init+366-leap; 
       leap = (((yr[nt]+1)%4+3)/4)
       last_val = last_val+366-leap
   return data_store        

# QUIVER ARROW FREQUENCY
afreq = 1
fSv = 1000000
Y2 = np.ones(len(mer_plt_time))

fig,ax = plt.subplots(1, figsize = (16,5))

Q = ax.quiver(mer_plt_time[::afreq],Y2[::afreq],Uwind[::afreq],Vwind[::afreq],width=0.001, scale=90)
qk = ax.quiverkey(Q,.05,.84,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})

ax.plot(ct_plt_time,pos_transport/fSv,'-r')
ax.plot(ct_plt_time,-1*neg_transport/fSv,'-b')
ax.set_xlim(mer_plt_time[0],ct_plt_time[-1])
ax.set_ylim(0,1.2)

ax.xaxis_date()

plt.title("Eastward (red) and Westward (blue) Transports in the upper 50m (Sv)")


# climatology and all yrs

dx = np.arange(0,365)
east_transp = clim_vals(pos_transport,1979,2007,ct_plt_time)
west_transp = clim_vals(neg_transport,1979,2007,ct_plt_time)

# Add winds - problem: starts

uwind = clim_vals(Uwind,1979,2007,mer_plt_time)
vwind = clim_vals(Vwind,1979,2007,mer_plt_time)

fig,ax = plt.subplots(1, figsize = (16,5))

for nt in range(east_transp.shape[0]):
    ax.plot(dx,east_transp[nt,:]/fSv,color= [1,.8,.8])
    ax.plot(dx,-1*west_transp[nt,:]/fSv,color=[.8,.8,1])

avg = np.nanmean(east_transp,axis=0)/fSv
ste = np.nanstd(east_transp/fSv,axis=0)/np.sqrt(east_transp.shape[0])

ax.plot(dx,avg,'-r',linewidth=3)
ax.plot(dx,avg + ste,'-r')
ax.plot(dx, avg - ste,'-r')

avg = np.nanmean(west_transp,axis=0)/fSv
ste = np.nanstd(west_transp/fSv,axis=0)/np.sqrt(west_transp.shape[0])

ax.plot(dx,-1*avg,'-b',linewidth=3)
ax.plot(dx,-1*avg + ste,'-b')
ax.plot(dx, -1*avg - ste,'-b')

Y2 = 0.9*np.ones(len(dx))

avg_v = np.mean(vwind,axis=0)
avg_u = np.mean(uwind,axis=0)


v_anom = vwind-avg_v
u_anom = uwind-avg_u

Q = ax.quiver(dx[::afreq],Y2[::afreq],avg_v[::afreq],avg_u[::afreq],width=0.001, scale=90)
qk = ax.quiverkey(Q,.05,.74,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})

plt.xticks([0,59,120,181,243,304],['Jan','Mar','May','Jul','Sept','Nov'])
plt.xlim(dx[0],dx[-1])
ax.set_ylim(0,1)
plt.title("Eastward (red) and Westward (blue) Transports in the upper 50m (Sv)")

fig,ax = plt.subplots(1, figsize = (16,5))

Uwind = u_anom.flatten()
Vwind = v_anom.flatten()

# remove leaps
init_ref = pltd.date2num(dt.datetime(1979,1,1,0,0))
last_ref = pltd.date2num(dt.datetime(2008,1,1,0,0))

itime = np.where((mer_plt_time>=init_ref) & (mer_plt_time<=last_ref))
mer_time2 = mer_plt_time[itime[0]]


for nt in range(1980,2007,4):
    i_val = np.in1d(mer_time2,pltd.date2num(dt.datetime(nt,2,29,12,0,0)),assume_unique=True,invert=True)
    mer_time2 = mer_time2[i_val]

Y2 = np.ones(len(mer_time2))

#Q = ax.quiver(Uwind[::afreq],Vwind[::afreq],width=0.001, scale=90)
Q = ax.quiver(mer_time2[::afreq],Y2[::afreq],Uwind[::afreq],Vwind[::afreq],width=0.001, scale=90)
qk = ax.quiverkey(Q,.05,.84,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})
ax.xaxis_date()
plt.show()




