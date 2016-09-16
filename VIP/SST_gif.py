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
plot_titles = ['CORTAD','CTROMS (5km)','VIP (500m)'] 
title_vals = ['CORTAD','CTROMS','VIP']
color_vals = ['-k','r','-b']

pr = [1,0,0]
pc = [0,0,1]

my_data=genfromtxt('soi_data.txt', usecols=np.arange(0,2))
my_data = my_data[542:589+1,:]
soi = my_data[:,1]
la_soi = np.ma.masked_where(soi<0,soi)
el_soi = np.ma.masked_where(soi>0,soi)

time_store = np.zeros((3,208))
sst_store = np.zeros((3,208))

place_time = convert_time(np.arange(35071.25,36527,7))
# INITIAL FIGURE
fig, ax = plt.subplots(2,2, figsize=(19,12))
fig.tight_layout(pad=3)
plt.subplots_adjust(right=0.90)
cbar_ax = fig.add_axes([0.93,0.12,0.02,.75])
lines = []
for nmod in range(len(models)):
    file = file_dir + models[nmod]
    fid = nc.Dataset(file)
    time = fid.variables[time_vals[nmod]][:]
    if nmod == 0:
       ind = np.arange(2,210)
       sst=get_cortad_vals(fid,ind)
       time = time[2:-1]
    else:
       sst = np.squeeze(fid.variables[sst_vals[nmod]][:])
    plot_times = convert_time(time) 
    pcol = ax[pr[nmod],pc[nmod]].pcolor(np.squeeze(sst[0,:,:]),vmin=23,vmax=32)
    ax[pr[nmod],pc[nmod]].set_title(plot_titles[nmod])
    ax[pr[nmod],pc[nmod]].set_axis_bgcolor((0.2,0.2,0.2))
    ax[pr[nmod],pc[nmod]].set_ylim(0,sst.shape[1])
    ax[pr[nmod],pc[nmod]].set_xlim(0,sst.shape[2])   

    domain_avg = np.mean(np.mean(sst,axis=2), axis=1) 
    lines.append(ax[1,1].plot(plot_times, domain_avg, color_vals[nmod],alpha=0.4,label=title_vals[nmod]))
    time_store[nmod,:] = plot_times
    sst_store[nmod,:] = domain_avg
    fid.close()

# LINE PLOT STUFF
ax[1,1].set_yticks([26,28,30])
ax[1,1].set_ylim(18,38)
ax[1,1].set_ylabel('oC')
ax[1,1].yaxis.set_label_coords(-.05,.5)

# SOI BAR PLOT
ax2 = ax[1,1].twinx()
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

ax2.text(place_time[93],.4,r'$El Ni\~{n}o$', color ='r', horizontalalignment='center')
ax2.text(place_time[155],-1,r'$La Ni\~{n}a$', color ='b', horizontalalignment='center')
ax2.set_yticks([-2,0,2])
ax2.set_ylim(-3,14)
ax2.set_ylabel('SOI',rotation=270)
ax2.yaxis.set_label_coords(1.05,.19)

ax[1,1].xaxis_date()
ax[1,1].xaxis.set_ticks_position('bottom')
ax[1,1].set_xlim(pltd.date2num(dt.datetime(1996,1,1)),pltd.date2num(dt.datetime(1999,12,31)))
handles, labels = ax[1,1].get_legend_handles_labels()
ax[1,1].legend(handles, labels, bbox_to_anchor=(0.98, 0.75) , ncol=3)
ax[1,1].text(place_time[3],37, 'MERRA VIP Winds', horizontalalignment='left')
ax[1,1].text(place_time[3],31.5, 'Domain-Averaged SST', horizontalalignment='left')

# WIND PLOT
wind_file = 'MERRA_VIP_weekly_1996-1999.nc'
fidw = nc.Dataset(wind_file)
timew = fidw.variables['time'][:]
datew = convert_time(timew)
u = fidw.variables['Uwind'][:]
v = fidw.variables['Vwind'][:]
avg_u = np.mean(np.mean(u,axis=2),axis=1)
avg_v = np.mean(np.mean(v,axis=2),axis=1)

Y2 = np.ones(len(avg_u))*35
afreq = 1

Q = ax[1,1].quiver(datew[::afreq],Y2[::afreq],avg_u[::afreq],avg_v[::afreq],width=0.001, scale=45)
qk = ax[1,1].quiverkey(Q,.8,.94,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':14})
fidw.close()
# FINAL DETAILS
cbar = plt.colorbar(pcol,cax=cbar_ax)
cbar_ax.text(-.3,1.04,'SST (oC)',multialignment='center')
#plt.show()

# MAKE THE ANIMATION
def get_sst(i,nmod):
     file = file_dir + models[nmod]
     fid = nc.Dataset(file)
     if nmod == 0:
       sst=get_cortad_vals(fid,i+2)
     else:
       sst = np.squeeze(fid.variables[sst_vals[nmod]][i,:])
     return sst
     fid.close()
gif_name = ['fifth.gif']
init = 160
endval = 208
for ntimes in range(2):
 ims = []
 n=0
 for add in range(init,endval):
    print add
    x_time = place_time[add]
    
    soi_val = check_dates-x_time
    soi_val2 = soi_val[abs(soi_val)<16]
    soi_time = soi_date[abs(soi_val)<16]
    el_val = el_soi[abs(soi_val)<16]
    la_val =  la_soi[abs(soi_val)<16]
    bar_width = bar_widths2[abs(soi_val)<16]

    if len(soi_val) > 1:
       if soi_val2[0]<bar_width[0]/2.:
          ind_val = 0
       else: 
          ind_val = 1
       soi_time = soi_time[ind_val]
       el_val = el_val[ind_val]
       la_val = la_val[ind_val]
       bar_width = bar_width[ind_val]
               
    sst1 = get_sst(add,0)
    sst2 = get_sst(add,1)
    sst3 = get_sst(add,2)
    lin2_ims, = ax[1,1].plot([time_store[0,add]],[sst_store[0,add]],'ko',markersize=5)
    lin3_ims, = ax[1,1].plot([time_store[1,add]],[sst_store[1,add]],'ro',markersize=5)
    lin4_ims, = ax[1,1].plot([time_store[2,add]],[sst_store[2,add]],'bo',markersize=5)
    bar_val1, = ax2.bar(soi_time,el_val,bar_width,color='r')
    bar_val2, = ax2.bar(soi_time,la_val,bar_width,color='b')
    lin_ims, = ax[1,1].plot([x_time,x_time],[10,40],'--g',linewidth=2)
    ims.append((ax[1,0].pcolor(sst1,vmin=23,vmax=32),\
                ax[0,0].pcolor(sst2,vmin=23,vmax=32),\
                ax[0,1].pcolor(sst3,vmin=23,vmax=32),\
                lin2_ims,lin3_ims,lin4_ims,\
                bar_val1, bar_val2,\
                lin_ims,))
 
 im_ani = animation.ArtistAnimation(fig, ims, blit=True)
 im_ani.save(gif_name[ntimes], writer = 'imagemagick',fps=10)
 init+=40; endval+=40
