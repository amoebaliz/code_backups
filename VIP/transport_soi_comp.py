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
tport_file = '/t3/workdir/liz/MODELS/VIP/PostProc/vip_transports/vip_u_ALLm_transport_col280_rows145-178.nc'
models = ['VIP_weekly_sst_1996-1999.nc']
sst_vals = ['temp']
time_vals = ['ocean_time']
plot_titles = ['VIP (500m)']
title_vals = ['VIP','Tran_temp']
color_vals = ['-b','-g']
e_lim = [60,200]
x_lim = [180,320]
#my_data=genfromtxt('soi_data.txt', usecols=np.arange(0,2))
#my_data = my_data[542:589+1,:]
#soi = my_data[:,1]
#la_soi = np.ma.masked_where(soi<0,soi)
#el_soi = np.ma.masked_where(soi>0,soi)

place_time = convert_time(np.arange(35071.25,36527,7))
# INITIAL FIGURE
fig, ax = plt.subplots(1, figsize=(19,12))
lines = []
for nmod in range(len(models)):
    file = '/t3/workdir/liz/MODELS/VIP/PostProc/' + models[nmod]
    fid = nc.Dataset(file)
    time = fid.variables[time_vals[nmod]][:]
    sst = np.squeeze(fid.variables[sst_vals[nmod]][:])

    cortad_fid = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/CORTAD_VIP_sst_1996-1999.nc')
    cortad_sst = cortad_fid.variables['sst'][:]
    cortad_sst = np.ma.masked_where(cortad_sst>100,cortad_sst)
    cortad_domain_avg = np.mean(np.mean(cortad_sst, axis=2),axis=1)
    cortad_trans_avg = np.mean(np.mean(cortad_sst[:,e_lim[0]:e_lim[1]+1,x_lim[0]:x_lim[1]+1],axis=2), axis=1)
#    cortad_trans_avg = np.mean(np.mean(cortad_sst[:,75:207,176:331],axis=2),axis=1)
    cortad_time = cortad_fid.variables['time'][:]
    ref = dt.datetime(1900,1,1,0,0)
    cortad_date_vals = np.zeros(len(cortad_time))
    for nt in range(len(cortad_time)):
        day_time = ref + dt.timedelta(days=np.float(cortad_time[nt]))
        cortad_date_vals[nt] = pltd.date2num(day_time)
    plot_times = convert_time(time)
    trans_sst = np.mean(np.mean(sst[:,e_lim[0]:e_lim[1]+1,x_lim[0]:x_lim[1]+1],axis=2),axis=1)
    avg_sst = np.mean(np.mean(sst,axis=2),axis=1)
#    ax.plot(plot_times, domain_avg, color_vals[nmod],alpha=0.4,label=title_vals[nmod])
    ax.plot(plot_times, trans_sst , color_vals[1],alpha=0.4,label=title_vals[1])
    ax.plot(plot_times, avg_sst,'-b') 
    ax.plot(cortad_date_vals,cortad_domain_avg,'-r')
    ax.plot(cortad_date_vals,cortad_trans_avg,'-k')

#    lines.append(ax.plot(plot_times, domain_avg, color_vals[nmod],alpha=0.4,label=title_vals[nmod]),ax.plot(plot_times, trans_sst , color_vals[1],alpha=0.4,label=title_vals[1]))
#    sst_store[nmod,:] = domain_avg
    fid.close()

# SST PLOT
ax.set_yticks([26,28,30])
ax.set_ylim(22,47)
ax.set_ylabel('oC')
ax.yaxis.set_label_coords(-.05,.43)

# TRANSPORT PLOT
fidT = nc.Dataset(tport_file)
time = fidT.variables['ocean_time'][:]
plot_time = convert_time(time)
pos_trans = fidT.variables['pos_transport'][:]
neg_trans = fidT.variables['neg_transport'][:]
ax3 = ax.twinx()
fSv = 1000000

ax3.plot(plot_time,pos_trans/fSv)
ax3.plot(plot_time,neg_trans/fSv)
#ax3.plot(plot_time,(pos_trans - neg_trans)/fSv,'-k')

# total_trans = pos_trans+neg_trans
# ax3.plot(plot_time,total_trans/fSv)
ax3.set_yticks([-.5,0,.5,1.0])
ax3.set_ylim(-2.5,2.3)
ax3.plot([plot_time[0],plot_time[-1]],[0,0],'-k')
ax3.set_ylabel('sV',rotation=270)
ax3.yaxis.set_label_coords(1.05,.7)

# SOI BAR PLOT
#ax2 = ax.twinx()
#bar_widths = [31,28,31,30,31,30,31,31,30,31,30,31]
#n = 0
#soi_date = np.zeros(4*12)
#bar_widths2 = np.zeros(4*12)
#check_dates = np.zeros(4*12)
#for nyr in range(1996,2000):
#    if nyr%4 == 0:
#       bar_widths[1] = 29
#    else:
#       bar_widths[1] = 28
#
#    for nmon in range(12):
#        date_val = dt.datetime(nyr,nmon+1,1)
#        check_dates[n] = pltd.date2num(date_val + dt.timedelta(days = bar_widths[nmon]/2. -1))
#        soi_date[n] = pltd.date2num(date_val)
#        bar_widths2[n] = bar_widths[nmon]
#        ax2.bar(soi_date[n],el_soi[n],bar_widths[nmon],color='r',alpha=0.4)
#        ax2.bar(soi_date[n],la_soi[n],bar_widths[nmon],color='b',alpha=0.4)

#        n+=1

#ax2.text(place_time[93],.4,r'$El Ni\~{n}o$', color ='r', horizontalalignment='center')
#ax2.text(place_time[155],-1,r'$La Ni\~{n}a$', color ='b', horizontalalignment='center')
#ax2.set_yticks([-2,0,2])
#ax2.set_ylim(-3,14)
#ax2.set_ylabel('SOI',rotation=270)
#ax2.yaxis.set_label_coords(1.05,.19)

ax.xaxis_date()
ax.xaxis.set_ticks_position('bottom')
ax.set_xlim(pltd.date2num(dt.datetime(1996,1,1)),pltd.date2num(dt.datetime(1999,12,31)))
handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, bbox_to_anchor=(0.7, 0.45) , ncol=1)
ax.text(place_time[3],38, 'MERRA VIP Winds', horizontalalignment='left')
ax.text(place_time[3],35.2, 'VIP Transport', horizontalalignment='left')
ax.text(place_time[3],30.5, 'Domain-Averaged SST', horizontalalignment='left')

# WIND PLOT
wind_file = '/t3/workdir/liz/external_data/MERRA/MERRA_VIP_wind.nc'
fidw = nc.Dataset(wind_file)
timew = fidw.variables['time'][6209:7670]
datew = np.zeros(len(timew))
for nt in range(len(timew)):
    day_time = ref + dt.timedelta(days=np.float(timew[nt]))
    datew[nt] = pltd.date2num(day_time)


u = fidw.variables['Uwind'][6209:7670,:,:]
v = fidw.variables['Vwind'][6209:7670,:,:]
avg_u = np.mean(np.mean(u,axis=2),axis=1)
avg_v = np.mean(np.mean(v,axis=2),axis=1)

Y2 = np.ones(len(avg_u))*43
afreq = 1

Q = ax.quiver(datew[::afreq],Y2[::afreq],avg_u[::afreq],avg_v[::afreq],width=0.001, scale=45)
qk = ax.quiverkey(Q,.8,.94,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})
# FINAL DETAILS
plt.show()
