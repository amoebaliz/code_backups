import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

def get_monthly_max(file):
    fid = nc.Dataset(file)
    temp = fid.variables['temp'][:]
    max_temp = np.squeeze(np.max(temp,axis=0))
    max_temp = np.ma.masked_where(max_temp>1000,max_temp)
    return max_temp    
   
    
fil_dir = '/t3/workdir/liz/MODELS/VIP/PostProc/monthly_temps/'

yrs = np.arange(1996,1999+1)
n=0
fid0 = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/monthly_temps/monthly_temps_1996.nc')
temp = fid0.variables['temp'][:]
max_vals = np.ma.zeros((len(yrs),temp.shape[1],temp.shape[2],temp.shape[3]))
print max_vals.shape
fid0.close()
for nyr in yrs:
    print nyr
    file = fil_dir + 'monthly_temps_' + str(nyr) + '.nc'
    max_vals[n,:,:,:] = get_monthly_max(file)
    n+=1

plt.figure()
plt.pcolor(np.squeeze(max_vals[0,49,:,:]))
plt.colorbar()
plt.show()
# Create/Save NetCDF file
ncid = nc.Dataset('/t3/workdir/liz/MODELS/CTROMS/VIP/monthly_temps/monthly_max_VIP_temps.nc','a')
ncid.createDimension('year', 20)

years = ncid.createVariable('year', 'f8',('year',))
years.units = 'year'
years[:] = yrs

max_temp = ncid.createVariable('temp', 'f8',('year','s_rho', 'eta_rho', 'xi_rho',))
max_temp.units = 'Celsius'
max_temp.long_name = 'time-averaged maximum monthly potential temperature'
max_temp.coordinates = 'lon_rho lat_rho s_rho ocean_time'
max_temp.field = 'temperature, scalar, series'
max_temp.time = 'year'
max_temp[:] = max_vals       
ncid.close()
