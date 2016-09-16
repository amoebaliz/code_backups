import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# file = 'CORTAD_VIP_sst_1996-1999.nc'
# temp_var = 'sst'
# mask_val = 1000

# file = 'CTROMS_VIP_sst_1996-1999.nc'
# temp_var = 'temp'

# file = 'VIP_sst_1996-1999.nc'
# temp_var = 'temp'

#-- old 2m sst --#

# VIP
# file = '/t3/workdir/liz/scripts/VIP_analyses/layer_sst/avg_2m_sst_96-99.nc'
# temp_var = 'layer_temp'

# CTROMS
#

#-- new VIP 2m sst --#




fid = nc.Dataset(file)
sst = fid.variables[temp_var][:].squeeze()
#sst = np.ma.masked_where(sst < mask_val, sst)
sst = np.ma.masked_where(np.isnan(sst), sst)

if np.ndim(sst)>2:
   sst= np.mean(sst,axis=0)

plt.figure()

plt.pcolor(sst)

plt.ylim(0,sst.shape[0])
plt.xlim(0,sst.shape[1])
plt.clim(27.5,29.5)
plt.colorbar()
plt.show()


