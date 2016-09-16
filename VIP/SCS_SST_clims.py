import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dt
from scipy.interpolate import interp1d

def index_weeks(yr):
    print yr
    tdelt_1 = (dt.datetime(yr,1,1,0,0)-ref).days
    tdelt_2 = (dt.datetime(yr+1,1,1,0,0)-ref).days
    itime = np.where((tdelt_1 <= cor_time) & (cor_time < tdelt_2))
    return itime[0]


# -@moeba---------------------------------- #

ncfile = '/data/external/P1/Data/CORTAD/Version4/remapped_CORTADv4_FilledSST.nc'
fid = nc.Dataset(ncfile)
cor_time = fid.variables['time'][:]
ref = dt.datetime(1900,1,1,0,0)
yrs = np.arange(1982,2009+1)
xyr = 1998 

# BOX BOUNDS
j0 = [460]
jN = [570]
i0 = [270]
iN = [320]

xnew = range(1,365+1)

for nreg in range(len(j0)):
    plt.figure()
    # STORAGE VAR
    dtim = yrs[-1]-yrs[0]
    dj = jN[nreg]-j0[nreg]+1
    di = iN[nreg]-i0[nreg]+1
    y_store = np.zeros((dtim*dj*di,365))
    n=0     
    for nyr in yrs:
        itime = index_weeks(nyr)
        an_temp = np.squeeze(fid.variables['sst'][itime[0]-1:itime[-1]+2,j0[nreg]:jN[nreg]+1,i0[nreg]:iN[nreg]+1])
        x_vals = cor_time[itime[0]-1:itime[-1]+2]-(dt.datetime(nyr,1,1,0,0)-ref).days

        if nyr == xyr:
           y_vals =  np.mean(an_temp,axis = (2,1))
           err_vals = np.std(an_temp,axis = (2,1))
           plt.errorbar(x_vals, y_vals, yerr = err_vals, fmt='ok', ecolor=[.5,.5,.5],elinewidth=3,capsize=0,zorder=30) 

        else:
           # Iterate/ interpolate over each spatial point 
           for jval in range(dj):
               for ival in range(di):
                   f = interp1d(x_vals,np.squeeze(an_temp[:,jval,ival]))
                   y_store[n,:] = f(xnew)
                   n+=1

    plt.errorbar(xnew,np.mean(y_store,axis=0),yerr=np.std(y_store,axis=0),fmt='-g',ecolor=[.8,.85,.75],elinewidth=2,capsize=0)

    title_str = 'eta = ' + str(j0[nreg]) + ':' + str(jN[nreg]) + ', xi = ' + str(i0[nreg]) + ':' + str(iN[nreg])
    plt.title(title_str)

plt.show()


        #y_store = np.zeros((yrs[-1]-yrs[0],365))
        #xnew = range(1,365+1)
        #n=0

        #an_temp = np.squeeze(fid.variables['sst'][itime[0]-1:itime[-1]+2,eta_val,xi_val])
        #x_vals = cor_time[itime[0]-1:itime[-1]+2]-(dt.datetime(nyr,1,1,0,0)-ref).days
        #f = interp1d(x_vals,an_temp)

        #   plt.plot(x_vals, an_temp)

