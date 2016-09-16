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

# OBJECTIVE: Map of the coral triangle with maximum DHW during 1998
# IDEAS: 1) Alternate Calculation of DHW
#        2) Inset map of the Philippines
#        3) Convert time to dates for the 'WHEN' plot

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(days=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

def polygon_patch(mapid,axs):
    mapid.drawcoastlines(linewidth=0)
    polys = []
    for polygon in mapid.landpolygons:
        polys.append(polygon.get_coords())

    lc = PolyCollection(polys, edgecolor='black',
         facecolor='0.8', closed=False)
    axs.add_collection(lc)

def plt_dhw(file):
    fid = nc.Dataset(fid)
 

##############################################
VIP = pyroms.grid.get_ROMS_grid('VIP')
deg_off = 0.02
# FILE DETAILS
dhw_files = ('/data/external/P1/Data/CORTAD/Version4/cortadv4_TSA_DHW_coral.nc',\
             '/t3/workdir/liz/MODELS/VIP/PostProc/DHW/Max_VIP_DHW_CoRTclim_1998.nc')
#             '/t3/workdir/liz/MODELS/VIP/PostProc/DHW/Max_VIP_DHW_CTclim_1998.nc')

cmap_file = '/home/frederic/python/cmap/dhw_noaa.cmap'
cmap = np.loadtxt(cmap_file)
lev_seven = np.array((159,212,0))
cmap2 = np.vstack((cmap[0:7,:],lev_seven,cmap[7:,:]))
cmap = cmap2/256.
dhw_noaa = colors.ListedColormap(cmap)

Ilats = [60,232]
Ilons = [168,335]

lats2 = VIP.hgrid.lat_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]
lons2 = VIP.hgrid.lon_rho[Ilats[0]:Ilats[1]+1,Ilons[0]:Ilons[1]+1]

cmin=0
cmax=16

# map domain
#lats2 = [12.4,15]
#lons2 = [120.,122.5]

for nt in range(len(dhw_files)):
    fid = nc.Dataset(dhw_files[nt])
    if nt < 1:
       # EXTRACT VARIABLES
       itime = range(844,895+1) # 1998
       time = fid.variables['time'][itime]
       plot_time = convert_time(time)
       lats = fid.variables['lat'][:]
       lons = fid.variables['lon'][:]

       DHW = fid.variables['dhw'][itime,:,:]
       DHW = np.ma.masked_where(DHW<0,DHW)
       DHW_max = np.max(DHW,axis=0)
       DHW_max_when = np.argmax(DHW,axis=0)
    elif nt < 2 :

       DHW_max = np.squeeze(fid.variables['max_DHW'][60:232,168:335])
       DHW_max_when = np.squeeze(fid.variables['year_day'][:])
       lats = fid.variables['lat_rho'][60:232,168:335]
       lons = fid.variables['lon_rho'][60:232,168:335]    
    else: 
       DHW_max = np.squeeze(fid.variables['max_DHW'][-1,:,:])
       DHW_max_when = np.squeeze(fid.variables['year_day'][:])
       lats = fid.variables['lat_rho'][:]
       lons = fid.variables['lon_rho'][:]


    # FIGURE DETAILS
    fig = plt.figure(figsize=(15,7))
    ax1 = fig.add_subplot(111)

    # Philippines Map
    m = Basemap(llcrnrlon=np.min(lons2)+.04,llcrnrlat=np.min(lats2)+.3,urcrnrlon=np.max(lons2)+deg_off,urcrnrlat=np.max(lats2)-.19,resolution='h',ax =ax1 )
    #m = Basemap(llcrnrlon=lons2[0],llcrnrlat=lats2[0],urcrnrlon=lons2[1],urcrnrlat=lats2[1],resolution='f')
    m.pcolormesh(lons,lats,DHW_max,vmin=cmin,vmax=cmax,cmap=dhw_noaa,latlon=True)
    polygon_patch(m,ax1)
    
#    m.drawmeridians([120.4,121.4])
#    m.drawparallels([13.4,13.9])


#    ax1.set_xlim(lons2[0], lons2[1])
#    ax1.set_ylim(lats2[0], lats2[1])

    #cb_ticks = np.linspace(cmin,cmax,5)
    #clb = m.colorbar(ticks=[])
#    for t in clb.ax.get_yticklabels():
#        t.set_fontsize(18)

    #plt.figure()
    #plt.pcolor(DHW_max_when,vmin=0,vmax=53)
    #plt.colorbar()
    fname = 'dhw_'+str(nt) +'.png'
    plt.savefig(fname)
plt.show()
