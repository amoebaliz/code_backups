import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# GRID FILE SPECIFICS
VIP = pyroms.grid.get_ROMS_grid('VIP')

dx = VIP.hgrid.dx
dy = VIP.hgrid.dy

# WIND STRESS FILE
#ncfile = '/t1/scratch/liz/tmpdir_VIP-LD.HCo13T/outputs/1996/VIP-LD.HCo13T_avg_1996-08-02T00:00:00.nc'
ncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/1996/VIP-LD.HCo13T_avg_1996-08-02T00:00:00.nc'

fid = nc.Dataset(ncfile)
taux = np.squeeze(fid.variables['sustr'][:])
tauy = np.squeeze(fid.variables['svstr'][:])
sst = np.squeeze(fid.variables['temp'][:,-1,:,:])
# Taux differenced along latitude/row axis
dTx = np.diff(taux,axis=0)
dy_u = (dy[:-1,:-1] + dy[:-1,1:] + dy[1:,:-1] + dy[1:,1:])/4.
# Tauy differenced along longitude/column axis
dTy = np.diff(tauy,axis=1)
dx_v = (dx[:-1,:-1] + dx[:-1,1:] + dx[1:,:-1] + dx[1:,1:])/4.
# Calculate wind stress curl
ws_curl = dTy/dx_v - dTx/dy_u

# Interpolate SST to PSI points
sst_2 = (sst[:-1,:-1] + sst[:-1,1:] + sst[1:,:-1] + sst[1:,1:])/4.
print sst_2.shape
print ws_curl.shape
fact = 1000000
c_range = 2

print np.max(ws_curl)
print np.min(ws_curl)

plt.figure()
plt.pcolor(fact*ws_curl,vmin = -1*c_range, vmax=c_range)
plt.colorbar()
plt.figure()
plt.pcolor(sst_2)
plt.show()
