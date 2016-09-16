import netCDF4 as nc
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd

fid = nc.Dataset('ctroms_u_transport_col395_monthly.nc','r')
pos_all = fid.variables['pos_transport'][:]
neg_all = -1*fid.variables['neg_transport'][:]
fid.close()



def plot_date(yrs):
    date = np.zeros(len(yrs)*12)
    nt = 0
    for nyr in yrs:
        for nmon in np.arange(1,12+1):
            # NOTE: the months are placed at the 1st of the month for simplicity
            day_time = dt.datetime(nyr,nmon,1,0,0)
            date[nt] = pltd.date2num(day_time)
            nt+=1
    return date

# 1982 and later
pos_sub = pos_all[264:]
neg_sub = neg_all[264:]

# RESHAPE and AVERAGE
p_sub = np.reshape(pos_sub,(-1,12)); p_clim = np.mean(p_sub,axis=0)
n_sub = np.reshape(neg_sub,(-1,12)); n_clim = np.mean(n_sub,axis=0)

fSv = 1000000
p_anom = (p_sub-p_clim)/fSv
n_anom = (n_sub-n_clim)/fSv

start = 1982
year1 = 1982; ind1 = year1-start
year2 = 2007; ind2 = year2-start
yrs =np.arange(year1,year2+1)
iyrs=yrs-start

p_anom_plot = np.ravel(p_anom[iyrs,:])
n_anom_plot = np.ravel(n_anom[iyrs,:])

plt_dates = plot_date(yrs)
print plt_dates.shape
print n_anom_plot.shape

fig,ax = plt.subplots(1)

ax.plot(p_clim/fSv,'-b')
ax.plot(n_clim/fSv,'-r')
ax.set_title('1982-2007 climatology')
plt.ylim(0,.22)

fig,ax = plt.subplots(1)

ax.plot(plt_dates,p_anom_plot,'-b')
ax.plot(plt_dates,n_anom_plot,'-r')
ax.plot([plt_dates[0],plt_dates[-1]],[0,0],'--k')

#ax.plot([pltd.date2num(dt.datetime(1996,12,1,0,0)),pltd.date2num(dt.datetime(1996,12,1,0,0))],[-.17,.17],'-g',linewidth=2)
#ax.plot([pltd.date2num(dt.datetime(1999,11,1,0,0)),pltd.date2num(dt.datetime(1999,11,1,0,0))],[-.17,.17],'-g',linewidth=2)
#ax.plot([pltd.date2num(dt.datetime(1996,12,1,0,0)),pltd.date2num(dt.datetime(1999,11,1,0,0))],[-.17,-.17],'-g',linewidth=2)
#ax.plot([pltd.date2num(dt.datetime(1996,12,1,0,0)),pltd.date2num(dt.datetime(1999,11,1,0,0))],[.17,.17],'-g',linewidth=2)




ax.set_title('1982-2007 transport anomalies')
ax.xaxis_date()
plt.ylim(-.2,.2)
plt.xlim(plt_dates[0],plt_dates[-1])
#fig,ax = plt.subplots(1,len(yrs),sharey=True)
#n = 0
#for nt in iyrs:
#    ax[n].plot(p_anom[nt,:]/fSv,'-b')
#    ax[n].plot(n_anom[nt,:]/fSv,'-r')
#    ax[n].plot([0,11],[0,0],'--k')
#    ax[n].set_title(str(yrs[n]))
#    ax[n].set_ylim(-.2,.2)
#    n+= 1
#plt.tight_layout(w_pad=0)
plt.show()
