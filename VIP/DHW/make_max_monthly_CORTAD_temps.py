import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dt
import pyroms

def index_weeks(yr):
    tdelt_1 = (dt.datetime(yr,1,1,0,0)-ref).days+3
    tdelt_2 = (dt.datetime(yr+1,1,1,0,0)-ref).days-3
    itime = np.where((tdelt_1 <= cor_time) & (cor_time < tdelt_2))
    return itime[0]

def save_ncfile(yr,max_vals):
    ncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/monthly_temps/monthly_max_sst_CORTAD_' + str(yr) + '.nc'
    ncid = nc.Dataset(ncfile,'a')
    ncid.createDimension('year', None)

    years = ncid.createVariable('year', 'f8',('year',))
    years.units = 'year'
    years[:] = yr

    max_temp = ncid.createVariable('max_temp', 'f8',('year', 'eta_rho', 'xi_rho',))
    max_temp.units = 'Celsius'
    max_temp.long_name = 'time-averaged maximum monthly potential temperature'
    max_temp.coordinates = 'lon_rho lat_rho year'
    max_temp.field = 'temperature, scalar, series'
    max_temp.time = 'year'
    max_temp[0,:] = max_vals
    ncid.close()
# ------------------------------------ #
GRID = pyroms.grid.get_ROMS_grid('VIP')
j,i = GRID.hgrid.mask_rho.shape

ncfile = '/t3/workdir/liz/MODELS/VIP/PostProc/CORTAD_to_VIP_remap/remapped_CORTADv4_FilledSST_to_VIP.nc'
fid = nc.Dataset(ncfile)
cor_time = fid.variables['time'][:]
ref = dt.datetime(1900,1,1,0,0)

yrs = np.arange(1982,2010+1)

for nyr in yrs:
    print nyr
    # storage variables 
    mon_max  = np.zeros((j,i))
    n=0

    itime = index_weeks(nyr)
    print itime[0],itime[-1],len(itime)
    for nwk in itime[:-3]:
        mon_temp = np.mean(fid.variables['sst'][nwk:nwk+4,:,:],axis=0)
        mon_temp = np.ma.masked_where(mon_temp > 1000, mon_temp)
        mon_max = np.where(mon_max<mon_temp,mon_temp,mon_max)
        mon_max = np.ma.masked_where(mon_max > 1000,mon_max)
    plt.figure()
    plt.pcolor(mon_max)
    plt.colorbar()
    plt.show()
    save_ncfile(nyr,mon_max)
