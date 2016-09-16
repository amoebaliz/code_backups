import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms

def get_clim(yr):

    cort_file = '/t3/workdir/liz/MODELS/VIP/PostProc/monthly_temps/monthly_max_sst_CORTAD_clim_82-10.nc'
    fid = nc.Dataset(cort_file)
    tn = yr-1982
    clim_temp = np.mean(fid.variables['max_temp'][:tn,:,:],axis=0)
    clim_temp = np.ma.masked_where(abs(clim_temp) > 1000, clim_temp)
    return clim_temp
    fid.close()
    
def get_temp(dir,n):

    ncfile = dir + 'VIP-LD.HCo13T_avg_' + str(n).zfill(3) + '.nc'
    print ncfile
    fid = nc.Dataset(ncfile,'r')
    temp = np.squeeze(fid.variables['temp'][:,-1,:,:])
    return temp
    fid.close()

def next_temp(dir,n):

    if n<nwkdy:
       temp = get_temp(dir,n+1)
    else:
       temp = get_temp(dir,n+1) - get_temp(dir,n-nwkdy+1)
    return temp

def save_ncfile(yr,max_vals,when_vals):
    ncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/DHW/Max_VIP_DHW_CoRTclim_' + str(yr) + '.nc'
    ncid = nc.Dataset(ncfile,'a')

    max_DHW = ncid.createVariable('max_DHW', 'f8',('eta_rho', 'xi_rho',))
    max_DHW.units = 'Degree heating weeks'
    max_DHW.long_name = 'Maximum accumulated thermal stress (DHW) signal'
    max_DHW.coordinates = 'lon_rho lat_rho'
    max_DHW.field = 'temperature, scalar, series'
    max_DHW[:] = max_vals

    max_when = ncid.createVariable('year_day', 'f8',('eta_rho', 'xi_rho',))
    max_when.units = 'Day Number'
    max_when.long_name = 'Day in year when accumulated maximum DHW stress'
    max_when.coordinates = 'lon_rho lat_rho'
    max_when.field = 'day_number, scalar, series'
    max_when[:] = when_vals
    ncid.close()


# ------------------------------------ #
GRID = pyroms.grid.get_ROMS_grid('VIP')
j,i = GRID.hgrid.mask_rho.shape
nwkdy = 12*7

for yr in np.arange(1996,1999+1):
    dir = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/filelinks/' + str(yr)+'/'

    # storage variables 
    sum_vals  = np.ma.zeros((j,i))   
    max_vals  = np.ma.zeros((j,i))
    when_vals = np.ma.zeros((j,i))

    if yr%4 == 0:
       nyrd = 366
    else:
       nyrd = 365

    for n in np.arange(nyrd):
    #for n in np.arange(83,90):
        sum_vals = sum_vals + next_temp(dir,n) 
        if n == nwkdy:
           max_vals = sum_vals
           when_vals.mask = sum_vals.mask
           when_vals[~when_vals.mask]=n
           print type(when_vals)
           print type(max_vals)
        elif n > nwkdy:
           when_vals = np.ma.where(max_vals<sum_vals,n,when_vals)
           max_vals  = np.ma.where(max_vals<sum_vals,sum_vals,max_vals)

    # subtract clim, divide by # days per week (7)
    max_DHW = (max_vals - nwkdy*get_clim(yr))/7.

    # Save NC files with max and when
    save_ncfile(yr,max_DHW,when_vals)

