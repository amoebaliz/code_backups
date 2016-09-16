import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats
from scipy import signal
from numpy import linalg

# REGRESSION FUNCTION DEFINED
def regmap_ts(field,ts):
    field[np.isnan(field)]=0
    
    regmap = np.zeros((field.shape[1],field.shape[2]))
    regmap_sig = np.zeros(regmap.shape)

    for y in range(field.shape[1]):
        for x in range(field.shape[2]):
            a = field[:,y,x]; b = ts
            slope, intercept, r_value, p_value, std_err = stats.linregress(b,a) 
            regmap[y,x] = slope

            if p_value < 0.01:
               regmap_sig[y,x] = slope
            else:
               regmap_sig[y,x] = 0
        print y 
    return regmap, regmap_sig

# FILE INFO
dir = '/t3/workdir/liz/scripts/MERRA/MERRA_SLP/'
#vip_ts_file = 'ctroms_u_transport_col395.nc'
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
#pos = ts_fid.variables['pos_transport'][7672:]
#neg = ts_fid.variables['neg_transport'][7672:]
#sum_trans = pos+neg
ct_zeta = np.squeeze(ts_fid.variables['zeta'][7672:])
ct_zeta2 = np.squeeze(ts_fid2.variables['zeta'][7672:])



ct_zeta_dif = ct_zeta - ct_zeta2
ts_fid.close()

# DETREND
SLP_test = signal.detrend(SLP,axis=0)
transport_ts = signal.detrend(ct_zeta_dif)

# NORMALIZE
SLP_test = (SLP_test - np.mean(SLP_test,axis=0))/np.std(SLP_test,axis=0)
transport_ts = (transport_ts - np.mean(transport_ts))/np.std(transport_ts)

# GENERATE 2D REGRESSION MATRIX
regmap, regmap_sig = regmap_ts(SLP_test,transport_ts)

# SAVING VARIABLES FOR FUTURE PLOTTING
np.save('regress_result',regmap)
np.save('sig_reg_result',regmap_sig)

# PLOTTING REGRESSION FIGURE
lons, lats = np.meshgrid(lons,lats)

# create new figure, axes instances.
fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])
# setup mercator map projection.
m = Basemap(llcrnrlon=95.,llcrnrlat=-12.,urcrnrlon=160.,urcrnrlat=32.,\
             resolution='l',projection='merc')

im1 = m.pcolormesh(lons,lats,regmap_sig,shading='flat',cmap=plt.cm.bwr,latlon=True,vmin=-0.5,vmax=0.5)
#m.plot([378,468])
m.drawcoastlines()
#m.fillcontinents()

# draw parallels
m.drawparallels(np.arange(-10,40,10),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
ax.set_title('Regression of SLP onto SSH dif (A-B), 1979-2007')

cb = m.colorbar(im1,"bottom", size="5%", pad="8%")
#cb.set_label('Regression Slope')

plt.show()
