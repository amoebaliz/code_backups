import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.animation as animation
import datetime as dt
from numpy import genfromtxt

def get_cortad_vals(fid,i):
    sst = fid.variables['sst'][i,:]
    sst2 = np.ma.masked_where(sst>300, sst)
    return sst2

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(days=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

# SET UP
file_dir = '/t3/workdir/liz/scripts/VIP_analyses/' 

models = ['CORTAD_VIP_sst_1996-1999.nc',\
          'CTROMS_VIP_weekly_sst_1996-1999.nc',\
          'VIP_weekly_sst_1996-1999.nc']

sst_vals = ['sst','temp','temp']
time_vals = ['time','ocean_time','ocean_time'] 
#plot_titles = ['CORTAD','CTROMS (5km)','VIP (500m)'] 

title_vals = ['CORTAD','CTROMS','VIP']
color_vals = ['-k','r','-b']

#pr = [1,0,0]
#pc = [0,0,1]

my_data=genfromtxt('soi_data.txt', usecols=np.arange(0,2))
my_data = my_data[542:589+1,:]
soi = my_data[:,1]
la_soi = np.ma.masked_where(soi<0,soi)
el_soi = np.ma.masked_where(soi>0,soi)

time_store = np.zeros((3,208))
sst_store = np.zeros((3,208))

# INITIAL FIGURE
fig, ax = plt.subplots(figsize=(12,8))
fig.tight_layout(pad=3)
#plt.subplots_adjust(right=0.90)
#cbar_ax = fig.add_axes([0.93,0.12,0.02,.75])
lines = []
#for nmod in range(len(models)):
for nmod in [2]:
    file = file_dir + models[nmod]
    fid = nc.Dataset(file)
    time = fid.variables[time_vals[nmod]][:]
    if nmod == 1:
       ind = np.arange(2,210)
       sst=get_cortad_vals(fid,ind)
       time = time[2:-1]
    else:
       sst = np.squeeze(fid.variables[sst_vals[nmod]][:])
    plot_times = convert_time(time) 
    domain_avg = np.mean(np.mean(sst,axis=2), axis=1) 
    lines.append(ax.plot(plot_times, domain_avg, color_vals[nmod],linewidth=2,label=title_vals[nmod]))
    lines.append(ax.plot(plot_times, sst[:,160,270], '-m',linewidth=2))
    time_store[nmod,:] = plot_times
    sst_store[nmod,:] = domain_avg
    fid.close()

# LINE PLOT STUFF
ax.set_yticks([26,28,30])
ax.set_ylim(20,32)
#ax.set_ylabel('oC')
ax.yaxis.set_label_coords(-.05,.5)

# SOI BAR PLOT
ax2 = ax.twinx()
bar_widths = [31,28,31,30,31,30,31,31,30,31,30,31]
n = 0 
soi_date = np.zeros(4*12)
bar_widths2 = np.zeros(4*12)
check_dates = np.zeros(4*12)
for nyr in range(1996,2000):
    if nyr%4 == 0:
       bar_widths[1] = 29
    else:
       bar_widths[1] = 28

    for nmon in range(12):
        date_val = dt.datetime(nyr,nmon+1,1)
        check_dates[n] = pltd.date2num(date_val + dt.timedelta(days = bar_widths[nmon]/2. -1))
        soi_date[n] = pltd.date2num(date_val)
        bar_widths2[n] = bar_widths[nmon]
        ax2.bar(soi_date[n],el_soi[n],bar_widths[nmon],color='r',alpha=0.4)
        ax2.bar(soi_date[n],la_soi[n],bar_widths[nmon],color='b',alpha=0.4)
        
        n+=1

#ax2.text(place_time[93],.4,r'$El Ni\~{n}o$', color ='r', horizontalalignment='center')
#ax2.text(place_time[155],-1,r'$La Ni\~{n}a$', color ='b', horizontalalignment='center')
ax2.set_yticks([-2,0,2])
ax2.set_ylim(-3,14)
#ax2.set_ylabel('SOI',rotation=270)
ax2.yaxis.set_label_coords(1.05,.19)

ax.xaxis_date()
ax.xaxis.set_ticks_position('bottom')
ax.set_xlim(pltd.date2num(dt.datetime(1996,1,1)),pltd.date2num(dt.datetime(1999,12,31)))
handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, bbox_to_anchor=(0.98, 0.75) , ncol=3)
#ax.text(place_time[3],31.5, 'Domain-Averaged SST', horizontalalignment='left')

plt.show()
