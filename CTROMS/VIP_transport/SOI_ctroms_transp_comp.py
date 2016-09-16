# Comparing with SOI values 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
import csv 

#### NOAA SOI RECORD
f = open('SOI_data.csv')
csv_f = csv.reader(f)
n = 0
soi_yr = []
soi_mon = []
soi = []
for row in csv_f:
   if n>1:
     soi_yr.append(int(row[0][0:4]))
     soi_mon.append(int(row[0][4:6]))
     soi.append(float(row[1]))
   n=n+1


print soi_yr

file = 'ctroms_u_transport_col395.nc'
fid = nc.Dataset(file,'r')

ref = dt.datetime(1900,1,1,0,0)
time = fid.variables['ocean_time'][:]
pos_transport = fid.variables['pos_transport']
neg_transport = fid.variables['neg_transport']

print pos_transport[0:5]
print neg_transport[0:5]

def monthly_means(ref,time_in,daily_data):
    dt_time = []
    for kt in range(len(time_in)):
        dt_time.append(ref + dt.timedelta(seconds=time_in[kt]))

    years = range(dt_time[0].year,dt_time[-1].year+1)
    monthly_means = np.zeros((len(years),12))    
    comp_time = pltd.date2num(dt_time)
    for year in years:

        for month in range(1,12+1):
            ref1 = pltd.date2num(dt.datetime(year,month,1,0,0))
            
            if month == 12:
                ref2 = pltd.date2num(dt.datetime(year+1,1,1,0,0))
            else:
                ref2 = pltd.date2num(dt.datetime(year,month+1,1,0,0))

            ivals = np.where((comp_time>=ref1)&(comp_time<=ref2))
            monthly_means[year-years[0],month-1] = np.mean(daily_data[ivals[0]])
            
    return monthly_means, years

pos_trans, years = monthly_means(ref,time,pos_transport)[:]
neg_trans = monthly_means(ref,time,neg_transport)[0]

print pos_trans.shape 
print len(years)







