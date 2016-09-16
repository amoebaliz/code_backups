import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats
from scipy import signal
from numpy import linalg

# REGRESSION FUNCTION DEFINED
def cormap_ts(field,ts):
    cormap = np.zeros((field.shape[1],field.shape[2]))

    for y in range(field.shape[1]):
        for x in range(field.shape[2]):
            a = field[:,y,x]; b = ts
            cor_mat = np.corrcoef(a,b) 
            cormap[y,x] = cor_mat[0,1]
        print y 
    return cormap

# FILE INFO
dir = '/t3/workdir/liz/scripts/MERRA/MERRA_SLP/'
vip_ts_file = '/t3/workdir/liz/scripts/ct_zeta_row468_col378.nc'
vip_ts_file2 = '/t3/workdir/liz/scripts/ct_zeta_row472_col402.nc'
files = dir + 'Pair_MERRA_daily_*.nc'
mfid = nc.MFDataset(files,'r')

# SLP INFO
lons = mfid.variables['lon'][145:235]
lats = mfid.variables['lat'][160:241]
mer_time = mfid.variables['time'][:-2194]
SLP = mfid.variables['Pair'][:-2194,160:241,145:235]
mfid.close()

# VIP SSH/TRANSPORT INFO 
#ilat = 
#ilon =

ts_fid = nc.Dataset(vip_ts_file,'r')
ts_fid2 = nc.Dataset(vip_ts_file2,'r')

ct_time = ts_fid.variables['ocean_time'][7672:]
ct_zeta = np.squeeze(ts_fid.variables['zeta'][7672:])
ct_zeta2 = np.squeeze(ts_fid2.variables['zeta'][7672:])
ct_zeta_dif = ct_zeta - ct_zeta2

ts_fid.close()


# GENERATE 2D REGRESSION MATRIX
cormap = cormap_ts(SLP,ct_zeta_dif)

# PLOTTING REGRESSION FIGURE
lons, lats = np.meshgrid(lons,lats)

# create new figure, axes instances.
fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])
# setup mercator map projection.
m = Basemap(llcrnrlon=95.,llcrnrlat=-12.,urcrnrlon=160.,urcrnrlat=32.,\
             resolution='l',projection='merc')

im1 = m.pcolormesh(lons,lats,cormap,shading='flat',cmap=plt.cm.bwr,latlon=True,vmin=-1,vmax=1)
#m.plot([378,468])
m.drawcoastlines()
#m.fillcontinents()

# draw parallels
m.drawparallels(np.arange(-10,40,10),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
ax.set_title('Correlation between SLP and SSH (A-B), 1979-2007')

cb = m.colorbar(im1,"bottom", size="5%", pad="8%")
#cb.set_label('Regression Slope')

plt.show()
