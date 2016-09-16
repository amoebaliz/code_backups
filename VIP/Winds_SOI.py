import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
#################################################################
# FILE INPUTS

# ----- WINDS ------ #
# Put in file names for winds; can be U/V separate or a single file

u_file = 'Uwinds.nc'
v_file = 'Vwinds.nc'

u_file = '/t3/workdir/liz/external_data/MERRA/Uwind_MERRA_daily_9699.nc'
v_file = '/t3/workdir/liz/external_data/MERRA/Vwind_MERRA_daily_9699.nc'


fidu = nc.Dataset(u_file)
fidv = nc.Dataset(v_file)
timew = fidu.variables['time'][:]

ref = dt.datetime(1900,1,1,0,0) # time reference
datew = np.zeros(len(timew))
for nt in range(len(timew)):
    # make sure to know if time is in days or seconds since reference point in below calculation
    day_time = ref + dt.timedelta(days=np.float(timew[nt]))
    datew[nt] = pltd.date2num(day_time)

# variables
u = fidu.variables['Uwind'][:]
v = fidv.variables['Vwind'][:]
avg_u = np.mean(np.mean(u,axis=2),axis=1)
avg_v = np.mean(np.mean(v,axis=2),axis=1)

# ----- INITIAL FIGURE ----- #
fig, ax = plt.subplots(1, figsize=(18,12))
# Edit date limits for x-axis
ax.xaxis.set_ticks_position('bottom')
ax.set_ylim(2,8) # arbitrary here - can be for specific variabile like temp or transport.
ax.set_yticks([])

# ---- WIND QUIVER PLOT ----------- #
ax2 = ax.twinx()
# Setting arbitrary Y-axis position for the winds relative to ax.set_ylim range above 
# form is: Y2 = np.ones(len(wind_data))*axis_position_scalar 
Y2 = np.ones(len(avg_u))*5
afreq = 1 # adjusts arrow frequency. example: 1 = every time point; 10 = every 10 time points 

Q = ax2.quiver(datew[::afreq],Y2[::afreq],avg_u[::afreq],avg_v[::afreq],width=0.001, scale=45)

# a quiverkey puts adds an arrow legend
qk = ax2.quiverkey(Q,.8,.94,1,r'$1 \frac{m}{s}$', labelpos='W')
ax2.set_yticks([])

# ---- OTHER PLOTS ---------------- #
# add other plots as needed
#ax3 = ax.twinx()
#ax3.set_yticks()
#ax3.set_ylim() # adjust for spacing relative to other plots
#ax3.plot()

# ----- FINAL DETAILS ------------- #
# edit your time bounds based on first and last date
ax.set_xlim(pltd.date2num(dt.datetime(1996,1,1)),pltd.date2num(dt.datetime(1999,12,31)))
# ax.set_xlim(pltd.date2num(dt.datetime(YYYY1,M1,D1)),pltd.date2num(dt.datetime(YYYYn,Mn,Dn)))
ax.xaxis_date()
#plt.savefig('INPUT_FIGURE_NAME.png')
plt.show()
