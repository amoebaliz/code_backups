import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.animation as animation
import datetime as dt
from numpy import genfromtxt

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(seconds=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

def plt_soi(ax):
    global data_stor
    my_data=genfromtxt('/t3/workdir/liz/external_data/NOAA_SOI/soi.txt', usecols=np.arange(1,13))
    my_data = my_data[45:48+1,:]
    soi = my_data.flatten()
    la_soi = np.ma.masked_where(soi<0,soi)
    el_soi = np.ma.masked_where(soi>0,soi)
    all_soi = []
    n = 0
    bar_widths = [31,28,31,30,31,30,31,31,30,31,30,31]
    soi_date = np.zeros(4*12)
    bar_widths2 = np.zeros(4*12)
    check_dates = np.zeros(4*12)
    ax.axhline(y=.5,color = 'cornflowerblue',ls='--')
    ax.axhline(y=-.5,color = 'salmon',ls='--')
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
            all_soi = np.hstack((all_soi,np.tile(soi[n],bar_widths[nmon])))

            ax.bar(soi_date[n],el_soi[n],bar_widths[nmon],color='salmon',zorder = 3)
            ax.bar(soi_date[n],la_soi[n],bar_widths[nmon],color='cornflowerblue',zorder = 3)
            n+=1
    ax.text(date_vals[370], -1, r'$El\ Ni\~{n}o$', color ='salmon', horizontalalignment='center')
    ax.text(date_vals[600], 1,r'$La\ Ni\~{n}a$', color ='cornflowerblue', horizontalalignment='center') 
    data_stor[nt,:] = all_soi

def plt_wnds(ax):
    global data_stor
    # WIND PLOT
    u_file = '/t3/workdir/liz/external_data/MERRA/Uwind_MERRA_daily_9699.nc'
    v_file = '/t3/workdir/liz/external_data/MERRA/Vwind_MERRA_daily_9699.nc'
    fidu = nc.Dataset(u_file)
    fidv = nc.Dataset(v_file)
    timew = fidu.variables['time'][:]
    datew = np.zeros(len(timew))
    for nt in range(len(timew)):
        day_time = ref + dt.timedelta(days=np.float(timew[nt]))
        datew[nt] = pltd.date2num(day_time)

    u = fidu.variables['Uwind'][:]
    v = fidv.variables['Vwind'][:]
    avg_u = np.mean(np.mean(u,axis=2),axis=1)
    avg_v = np.mean(np.mean(v,axis=2),axis=1)
    fil_val = np.random.rand(len(date_vals))/2.0
    data_stor[nplts-1,:] = fil_val  
    Y2 = np.ones(len(avg_u))*np.mean(fil_val)
    afreq = 1

    Q = ax.quiver(datew[::afreq],Y2[::afreq],avg_u[::afreq],avg_v[::afreq],width=0.001, scale=40)
    qk = ax.quiverkey(Q,.81,.98,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':18})
    ax.set_yticks([])
    ax.text(txpos, Y2[0]+1.4*np.std(fil_val), 'MERRA Winds', fontsize=16)
    data_stor[nplts-1,:] = fil_val 
    
def set_yaxes():
    # PADDING CALCULATIONS
    pltmin = np.mean(data_stor[0,:]) - stdfrac*np.std(data_stor[0,:])
    submax = np.mean(data_stor[0,:]) + stdfrac*(2*(nplts-1)+1)*np.std(data_stor[0,:])

    pltmax = np.mean(data_stor[nplts-1,:]) + stdfrac*np.std(data_stor[nplts-1,:])
    submin = np.mean(data_stor[nplts-1,:]) - stdfrac*(2*(nplts-1)+1)*np.std(data_stor[nplts-1,:])

    lpd = (1.0*pltmin-np.min(data_stor[0,:]))/(submax-pltmin) + pdad
    upd = (1.0*np.max(data_stor[nplts-1,:])-pltmax)/(pltmax-submin) + pdad

    # ITERATE OVER PLOTS    
    for nt in range(nplts):
        plt_rng = nplts*2*stdfrac*np.std(data_stor[nt,:])
        ymin = np.mean(data_stor[nt,:]) - stdfrac*(2*nt+1)*np.std(data_stor[nt,:])
        ymax = np.mean(data_stor[nt,:]) + stdfrac*(2*(nplts - (nt+1))+1)*np.std(data_stor[nt,:])

        # APPLY PADS
        ymin = ymin - lpd*plt_rng
        ymax = ymax + upd*plt_rng
        axs[nt].set_ylim(ymin,ymax)

        # LABEL POSITION CENTERED AT 0.0 TICK 
        if nt < nplts-1:
           ypos = (0.0 - ymin)/(ymax-ymin)
           axs[nt].set_yticks([-1*tcks[nt],0,tcks[nt]])
           axs[nt].set_ylabel(ylabels[nt], y = ypos, rotation = 90 + (nt%2)*180, fontsize=18, labelpad = 15)    
           axs[nt].tick_params(which='major', labelsize=18)
           if nt%2==1:
              axs[nt].yaxis.tick_right()
              axs[nt].yaxis.set_label_position("right")

           else:
              axs[nt].yaxis.tick_left()
              axs[nt].yaxis.set_label_position("left")
            
    #ax.figure.canvas.draw() 
# ------------------------------------------------------ #
# DATA DETAILS
dir = '/t3/workdir/liz/MODELS/VIP/PostProc/'
fils = ['vip_transports/vip_u_100m_transport_col278_rows144-178.nc',\
        'SSH_1996-1999.nc']#,\
#        'u_surf_278_96_99.nc']

vars = ['_transport','zeta','u']

# CONSTANTS
fSv = 1000000
ref = dt.datetime(1900,1,1,0,0)

# PLOT SPECS
months = pltd.MonthLocator(bymonth=7,bymonthday=1)
monthsFmt = pltd.DateFormatter('%B %Y')

nplts = len(fils)+2
colval = ['-b','-g','-k','-r']
ylabels = []
tcks = [2,.4,2,.5,0.01]
txpos = pltd.date2num(dt.datetime(1996,1,1,0,0))+30
xmin  = pltd.date2num(dt.datetime(1996,1,1,0,0))
xmax  = pltd.date2num(dt.datetime(1999,12,31,0,0))
ylabels = ['SOI','Transport (Sv)',r'$\Delta$ SSH (cm)','Velocity (m/s)','thing2 ()']
stdfrac = 2.5 # # OF STDEVs (*.5) DEVOTED TO PLOT SPACE
pdad = .01


fig, ax = plt.subplots(1, figsize=(12,14))
[fig.add_axes(ax.twinx()) for n in range(nplts-1)]
axs = fig.get_axes()


# PLOTTING and STORING
for nt in range(nplts-2):
    ncfile = dir + fils[nt]
    fid = nc.Dataset(ncfile)

    if nt == 0:
       # Get time for shared x-axis
       time = fid.variables['ocean_time'][:]
       date_vals = convert_time(time)

       # Initialize data storage matrix 
       data_stor = np.empty([nplts, len(date_vals)])

       # Plot SOI as first plot
       plt_soi(axs[nt])

       # TRANSPORT variable modifications
       pos_trans =  np.squeeze(fid.variables['pos_transport'][:])/fSv
       neg_trans =  np.squeeze(fid.variables['neg_transport'][:])/fSv

       # Plot transports as second plot       
       axs[nt+1].plot(date_vals,pos_trans,color='darkviolet')
       axs[nt+1].plot(date_vals,neg_trans,color='darkorange')
       
       var = (pos_trans+neg_trans)

# CONTINUE PLOTTING FOR SSH/ SURF VEL / WINDS
    else:
       var = np.squeeze(fid.variables[vars[nt]][:])
       if nt == 1:
          # SSH variable modification
          var = 100*(np.squeeze(np.mean(var[:,65:149,180],axis=1))-np.squeeze(np.mean(var[:,149:199,318],axis=1)))
       else:
          # Vertical variable modification
          var = np.mean(np.mean(var,axis=2),axis=1)

       axs[nt+1].plot(date_vals,var,colval[nt])

    axs[nt+1].axhline(y=0,color = 'k',zorder=0)
    data_stor[nt+1,:] = var

data_stor = data_stor.reshape(nplts,-1)

plt_wnds(axs[nplts-1])

set_yaxes()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)

# FINAL DETAILS
ax.xaxis_date()
ax.set_xlim(xmin,xmax)
ax.xaxis.set_ticks_position('bottom')

#plt.tight_layout(pad=2)
plt.savefig('drivers.png')
plt.show()
