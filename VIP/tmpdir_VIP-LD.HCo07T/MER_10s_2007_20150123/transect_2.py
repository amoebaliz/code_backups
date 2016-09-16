import netCDF4 as nc
import numpy as np
import pyroms as py
from matplotlib import pyplot as plt

data_file = 'JAN_current.nc'
fid = nc.Dataset(data_file,'r')

dstname = 'VIP'
dst_grd = py.grid.get_ROMS_grid(dstname)

u = fid.variables['u'[:]]
u_avg = np.mean(u[:,:,:,:], axis = 0)

slice,Z,X,Y = py.tools.islice(u_avg,300,dst_grd,'u')


plt.contourf(Y,Z,slice[::-1,:],100)
plt.show()


# import pyroms
# pyroms.<TAB>
