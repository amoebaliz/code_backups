import numpy as np
import netCDF4 as nc
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(days=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

def index_weeks(yr):
    tdelt_1 = (dt.datetime(yr,1,1,0,0)-ref).days
    tdelt_2 = (dt.datetime(yr+1,1,1,0,0)-ref).days
    itime = np.where((tdelt_1 <= cor_time) & (cor_time < tdelt_2))
    return itime[0]

def plt_sst_ts(ax):
    # STORAGE VAR
    dtim = yrs2[-1]-yrs2[0]
    dj = jN-j0+1
    di = iN-i0+1
    y_store = np.zeros((dtim*dj*di,365))
    n=0

    for nyr in yrs2:
        itime = index_weeks(nyr)

        an_temp = np.squeeze(sstfid.variables['sst'][itime[0]-1:itime[-1]+2,j0:jN+1,i0:iN+1])
        an_temp[an_temp<0]=np.NAN

        x_vals = cor_time[itime[0]-1:itime[-1]+2]-(dt.datetime(nyr,1,1,0,0)-ref).days
        if nyr == YOI:
           y_vals =  np.nanmean(an_temp,axis = (2,1))
           err_vals = np.nanstd(an_temp,axis = (2,1))
           xdate = convert_time(x_vals)
           ax.errorbar(xdate[1:-1], y_vals[1:-1], yerr = err_vals[1:-1], fmt='ok', ms = 2, ecolor=[.5,.5,.5],elinewidth=2,capsize=0,zorder=30)
        else:
           # Iterate/ interpolate over each spatial point 
           for jval in range(dj):
               for ival in range(di):
                   f = interp1d(x_vals,np.squeeze(an_temp[:,jval,ival]))
                   y_store[n,:] = f(xnew)
                   n+=1
    
    xdate=convert_time(xnew)
    ax.errorbar(xdate,np.nanmean(y_store,axis=0),yerr=np.nanstd(y_store,axis=0),fmt='-g',ecolor=[.8,.85,.75],elinewidth=2,capsize=0)

# ---------------------------------------------------------- #
months = pltd.MonthLocator(interval=2)
monthsFmt = pltd.DateFormatter('%b')

mer_dir = '/t3/workdir/liz/external_data/MERRA/MERRA_WINDS/'
sstfile = '/data/external/P1/Data/CORTAD/Version4/cortadv4_FilledSST_coral.nc'
sstfid = nc.Dataset(sstfile)
cor_time = sstfid.variables['time'][:]
ref = dt.datetime(1900,1,1,0,0)

# EXTRACT CORTAD VARIABLES
sstfid = nc.Dataset(sstfile)
cor_time = sstfid.variables['time'][:]
YOI = 1998
yrs2 = np.arange(1982,YOI+1)
xnew = range(1,365+1)

# VIP REGION BOUNDS
j0 = 860
jN = 875
i0 = 735
iN = 750

yrs = range(1979,2013+1)
U_store = np.zeros((len(yrs),365))
V_store = np.zeros((len(yrs),365))

for nt in range(len(yrs)):

    uncfile = mer_dir + 'Uwind_MERRA_daily_' + str(yrs[nt]) + '.nc'
    vncfile = mer_dir + 'Vwind_MERRA_daily_' + str(yrs[nt]) + '.nc'

    ufid = nc.Dataset(uncfile)
    vfid = nc.Dataset(vncfile)

    Uwind = np.mean(np.mean(ufid.variables['Uwind'][:],axis=2),axis=1)
    Vwind = np.mean(np.mean(vfid.variables['Vwind'][:],axis=2),axis=1) 

    U_store[nt,:] = Uwind[:365]
    V_store[nt,:] = Vwind[:365]

    if yrs[nt] == 1998:
       U_98 = Uwind
       V_98 = Vwind

# PLOTTING SPECIFICATIONS
X = range(365)
datew = convert_time(X)
Y1 = 3*np.ones(len(X))

xmin  = pltd.date2num(dt.datetime(1899,11,20,0,0))
xmax  = pltd.date2num(dt.datetime(1901,1,2,12,0))
ymin  = 24.8
ymax  = 36
y2min = -20.5
y2max = 7

avg_u = np.mean(U_store,axis=0)
avg_v = np.mean(V_store,axis=0)

afreq = 1

fig, ax = plt.subplots(1, figsize=(7,7))
ax2 = ax.twinx()

plt_sst_ts(ax)


Q1 = ax2.quiver(datew[::afreq],Y1[::afreq],avg_u[::afreq],avg_v[::afreq],width=0.0015, scale=20,clip_on=False)
qk1 = ax2.quiverkey(Q1,.9,.94,1,r'$1 \frac{m}{s}$', labelpos='W',fontproperties={'size':16})

UA_98 = U_98-avg_u
VA_98 = V_98-avg_v

# Just JJAS
UA_98  = UA_98[151:273]
VA_98  = VA_98[151:273]
datew2 = datew[151:273]
Y2 = -1.5*np.ones(len(UA_98))

#print np.mean(UA_98)
#print np.mean(VA_98)

Q2 = ax2.quiver(datew2[::afreq],Y2[::afreq],UA_98[::afreq],VA_98[::afreq],color='b',width=0.0015, scale=20,clip_on=False)
#qk2 = ax.quiverkey(Q2,.9,.14,1,r'$1 \frac{m}{s}$', labelpos='W',labelcolor='b',fontproperties={'size':18})

# AXES
ax.tick_params(axis='both', which='both', top='off',labelsize=16)

# X-AXIS 
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
ax.set_xlim(xmin,xmax)
ax.xaxis.set_ticks([pltd.date2num(dt.datetime(1900,1,1,0,0)),\
                   pltd.date2num(dt.datetime(1900,3,1,0,0)),\
                   pltd.date2num(dt.datetime(1900,5,1,0,0)),\
                   pltd.date2num(dt.datetime(1900,7,1,0,0)),\
                   pltd.date2num(dt.datetime(1900,9,1,0,0)),\
                   pltd.date2num(dt.datetime(1900,11,1,0,0))])
ax.xaxis.set_ticks_position('bottom')

# Y-AXIS
ax.set_ylim(ymin,ymax)
ax.yaxis.set_ticks([26,28,30])
ax.set_ylabel(u'Temperature (\N{DEGREE SIGN}C)', y= .29,labelpad = 20, fontsize=16)

ax2.set_ylim(y2min,y2max)
ax2.yaxis.set_ticks([])

# TEXT LABELS
txpos = pltd.date2num(dt.datetime(1900,1,1,0,0))-30
ax.text(txpos, 30.2, 'CoRTAD Climatology', color = 'g', fontsize=16)
ax.errorbar(txpos, 29.2, yerr = .5, fmt='ok', ms = 2, ecolor=[.5,.5,.5],elinewidth=2,capsize=0)
ax.text(txpos+5, 29.2, '1998', va = 'center', color = 'k', fontsize=16)

ax2.text(txpos, Y1[0]+2, 'MERRA Climatology', fontsize=16)
ax2.text(txpos, Y2[0], '1998 Anomaly', color = 'b', fontsize=16)

# FINAL DETAILS
plt.tight_layout(pad=2)
plt.savefig('wind.png')
plt.show()
