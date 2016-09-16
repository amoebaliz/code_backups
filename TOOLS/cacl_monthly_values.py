import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import netCDF4 as nc

# BASIC DETAILS
yfirst = 1960
ylast = 2007
infile = 'ctroms_zeta.nc'
var1_name = 'zeta'
#var2_name = 'neg_transport'
time_name = 'ocean_time'
outfile = 'ctroms_zeta_monthly.nc'

# VECTORS AND WHATNOT
yrs = range(yfirst,ylast)
fSv = 1000000
# FILE INFORMATION
fid = nc.Dataset(infile,'r')
daily_time = fid.variables[time_name]
daily_var1 = fid.variables[var1_name]
#daily_var2 = fid.variables[var2_name]
# WRITE NEW FILE

fid2 = nc.Dataset(outfile,'w')
fid2.createDimension(time_name,None)
time = fid2.createVariable(time_name,'f8',(time_name,))
time.units = daily_time.units

var1 = fid2.createVariable(var1_name,'f8',(time_name,))
var1.units = daily_var1.units
var1.long_name = daily_var1.long_name

#ADDITIONAL VARIABLES
#var2 = fid2.createVariable(var2_name,'f8',(time_name,))
#var2.units = daily_var2.units
#var2.long_name = daily_var2.long_name

daily_time = daily_time[732:]
daily_var1 = daily_var1[732:]
#daily_var2 = daily_var2[732:]

# TAKE DAILY MEANS AND AVERAGE OVER MONTHS
days_n_months = [31,28,31,30,31,30,31,31,30,31,30,31]
month_var1 = np.zeros(12*len(yrs),daily_var1.shape[1],daily_var1.shape[2])
#month_var2 = np.zeros(12*len(yrs))
month_time = np.zeros(12*len(yrs))

d_0 = 0
for nyr in yrs:
#for nyr in range(1960,1965):
    for nmon in range(12):
    # IF LEAP YEAR
        if nyr%4 == 0 and nmon == 1:
           ndays = 29
        else:
           ndays = days_n_months[nmon]

        # calculate the average value
        d_end = d_0 + ndays
     
        if nyr == ylast and nmon == 11:
        # just average to end of month
           month_var1[(nyr-yfirst)*12+nmon] = np.mean(daily_var1[d_0:,:,:],axis=0)
 #          month_var2[(nyr-yfirst)*12+nmon] = np.mean(daily_var2[d_0:])
           month_time[(nyr-yfirst)*12+nmon] = np.mean(daily_time[d_0:])
        else: 
           # calculate average from d_0:d_end
           month_var1[(nyr-yfirst)*12+nmon] = np.mean(daily_var1[d_0:d_end,:,:],axis=0)
 #          month_var2[(nyr-yfirst)*12+nmon] = np.mean(daily_var2[d_0:d_end])
           month_time[(nyr-yfirst)*12+nmon] = np.mean(daily_time[d_0:d_end])
        print month_var1[(nyr-yfirst)*12+nmon]/fSv   

        d_0 = d_end

time[:] = month_time
var1[:] = month_var1
#var2[:] = month_var2

fid2.close()

