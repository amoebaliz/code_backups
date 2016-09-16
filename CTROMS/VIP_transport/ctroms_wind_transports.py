# TIME SERIES OF MERRA WINDS ON AND CTROMS TRANSPORT THROUGH THE VIP
import netCDF4 as nc
import numpy as np

merra_dir = '/t1/scratch/forcings_sets/MERRA/drowned/'

min_lon = 120.5;  max_lon = 122.0
min_lat = 13;   max_lat = 14



def avg_wind(i,dir):

    uwind = merra_dir + 'drowned_MERRA_Uwind_3hours_' + str(i) + '.nc'
    vwind = merra_dir + 'drowned_MERRA_Vwind_3hours_' + str(i) + '.nc'
    
    Uwind = ufid.variables['Uwind'][:,Ilat[0],Ilon[0]]
    Vwind = vfid.variables['Vwind'][:,Ilat[0],Ilon[0]]


    var = var + 'wind'
    file = merra_dir + 'drowned_MERRA_' + var + 'wind_3hours_' + str(i) + '.nc'
    fid = nc.Dataset(file,'r')
    wind = fid.variables[var][:,Ilat[0],Ilon[0]]



    return uavg, vavg, time

for nt in range(1979,2013+1)
    if nt == 1979:
       # 

    u, time = avg_wind(nt)

    # open netcdf file in append mode
    # append u, 



ctroms_dir = ''

# PROCESS EACH YEAR OF MERRA WINDS: 
  # OPEN U and V, AVERAGE BY INCREMENTS OF 8, AND OVER VIP AREA
  # STORE in .nc file/ update it every year
  # Write a comment about the location


#PROCESS EVERY DAY IN THE CT:
  # OPEN U and V calculate transport at specific location
  # store transport and time data file. 
  # Write a comment about the location


# WRITE NCFILE DETAILS
# variables = ocean time, west transport, east transport
# global attribute

ct_nums = range(1,18262+1)

for nt in ct_nums:
    # get the pos and neg transport
    




