import pyroms
import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import shapefile as shp

# FIGURE OUT PERCENT COVERAGE IN CORTAD

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0)
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor='0.8', closed=False)
    axs.add_collection(lc)

 

##############################################
# FILE DETAILS
cortfile = ('/t3/workdir/liz/external_data/CORTAD/cortadv4_WeeklySST_coral.nc')
/data/external/P1/Data/CORTAD/Version4

fid = nc.Dataset(cortfile)

j0 = [860]
jN = [875]
i0 = [735]
iN = [750]

lats = fid.variables['lat'][j0:jN,i0:iN]
lons = fid.variables['lon'][j0:jN,i0:iN]

sst1 = fid.variables['sst'][itime1,j0:jN,i0:iN]
sst2 = fid.variables['sst'][itime2,j0:jN,i0:iN]

sst1 = np.ma.masked_where(sst1<0,sst1)
sst2 = np.ma.masked_where(sst2<0,sst2)

itime1 = range(843,896+1) # 1998
itime2 = range(864,883+1) # JJAS 1998

# BIG domain
Ilats = [60,232]
Ilons = [168,335]

# sub domain
lats2 = [12.4,15]
lons2 = [120.,122.5]

perct_stor1 = np.zeros((sst1.shape[1],sst1.shape[2]))

for j in range(sst1.shape[1]):
    for i in range(sst1.shape[2]):
        perct_stor1[j,i] = sst1[:,j,i].count()/(1.0*sst1.shape[0])        

print np.min(perct_stor1)



fname = 'percent_temps.png'


fig = plt.figure()
ax = fig.add_subplot(111)
m2 = Basemap(llcrnrlon=lons2[0],llcrnrlat=lats2[0],urcrnrlon=lons2[1],urcrnrlat=lats2[1],resolution='f')

m2.pcolor(lons,lats,perct_stor1)
polygon_patch(m2,ax)
plt.colorbar()

plt.figure()
plt.pcolor(perct_stor1)
plt.colorbar()
plt.savefig(fname)
plt.show()
