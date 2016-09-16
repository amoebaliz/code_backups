import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as pltd
import matplotlib.colors as colors
from scipy.interpolate import interp1d
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import shapefile as shp

# OBJECTIVE: Map of the coral triangle with maximum DHW during 1998
# IDEAS: 1) Alternate Calculation of DHW
#        2) Inset map of the Philippines
#        3) Convert time to dates for the 'WHEN' plot

def index_weeks(yr):
    tdelt_1 = (dt.datetime(yr,1,1,0,0)-ref).days
    tdelt_2 = (dt.datetime(yr+1,1,1,0,0)-ref).days
    itime = np.where((tdelt_1 <= cor_time) & (cor_time < tdelt_2))
    return itime[0]

def convert_time(time):
    ref = dt.datetime(YOI-1,12,31,0,0)
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

def plt_shp_coord(sf,m_fin):
    m0 = Basemap(projection = 'cea')
    x_0,y_0 = m0(0,0)
    for shape in sf.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        x = np.array(x)-x_0
        y = np.array(y)+y_0
        lons,lats = m0(x,y,inverse=True)
        m_fin.plot(lons,lats,'.1',lw=6)
        m_fin.plot(lons,lats,'paleturquoise',lw=1.5)

def plt_inset(ax):
    axins = zoomed_inset_axes(ax, 10, loc=1, borderpad=0)

    axins.set_xlim(lons2[0], lons2[1])
    axins.set_ylim(lats2[0], lats2[1])

    m2 = Basemap(llcrnrlon=lons2[0],llcrnrlat=lats2[0],urcrnrlon=lons2[1],urcrnrlat=lats2[1],resolution='f',ax=axins)
    m2.pcolormesh(lons,lats,DHW_max,vmin=0,vmax=16,cmap=dhw_noaa,latlon=True)
    polygon_patch(m2,axins)
    m2.drawmapboundary(linewidth=1.0)

    mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0",lw =1)
    return axins
    #gs1.tight_layout(fig, rect=[.03, .05, .8, .8])
    #pos = ax1.get_position()



def plt_sst_ts(nt,ax):
    # STORAGE VAR
    dtim = yrs[-1]-yrs[0]
    dj = jN[nt]-j0[nt]+1
    di = iN[nt]-i0[nt]+1
    y_store = np.zeros((dtim*dj*di,365))
    n=0

    #for nyr in yrs[:-1]:
    for nyr in yrs[:]:
        itime = index_weeks(nyr)

        an_temp = np.squeeze(sstfid.variables['sst'][itime[0]-1:itime[-1]+2,j0[nt]:jN[nt]+1,i0[nt]:iN[nt]+1])
        an_temp[an_temp<0]=np.NAN

        x_vals = cor_time[itime[0]-1:itime[-1]+2]-(dt.datetime(nyr,1,1,0,0)-ref).days
        if nyr == YOI:
        #if nyr == 0:
           y_vals =  np.nanmean(an_temp,axis = (2,1))
           err_vals = np.nanstd(an_temp,axis = (2,1))
           xdate = convert_time(x_vals)
           ax.errorbar(xdate[1:-1], y_vals[1:-1], yerr = err_vals[1:-1], fmt='ok', ms = 4, ecolor=[.5,.5,.5],elinewidth=2,capsize=0,zorder=30)

        else:
           # Iterate/ interpolate over each spatial point 
           for jval in range(dj):
               for ival in range(di):
                   f = interp1d(x_vals,np.squeeze(an_temp[:,jval,ival]))
                   y_store[n,:] = f(xnew)
                   n+=1

    xdate=convert_time(xnew)
    ax.errorbar(xdate,np.nanmean(y_store,axis=0),yerr=np.nanstd(y_store,axis=0),fmt='-g',ecolor=[.8,.85,.75],elinewidth=2,capsize=0)
    ax.set_xlim(xdate[0]-1,xdate[-1]+1)
    ax.xaxis.tick_bottom()
    ax.xaxis.set_major_locator(months)
    ax.set_xticklabels(labels)

    ax.yaxis.tick_right()
    ax.set_ylim(ymins[nt],ymaxs[nt])
    ax.set_yticks((yt1[nt],yt3[nt]))
    if nt == len(j0)-1:
       ax.set_xlabel('Time', fontsize=14, labelpad=10)
   
    yloc = (ymaxs[nt]-ymins[nt])/7.
    #ax.text(xdate[7],ymaxs[nt]-yloc,reg_labels[nt])
    
def plt_box(nt,ax):
 
    ax.plot((lons[j0[nt],i0[nt]],lons[j0[nt],iN[nt]]),(lats[j0[nt],i0[nt]],lats[j0[nt],iN[nt]]),'-w',lw=1.5)
    ax.plot((lons[j0[nt],i0[nt]],lons[jN[nt],i0[nt]]),(lats[j0[nt],i0[nt]],lats[jN[nt],i0[nt]]),'-w',lw=1.5)
    ax.plot((lons[j0[nt],iN[nt]],lons[jN[nt],iN[nt]]),(lats[j0[nt],iN[nt]],lats[jN[nt],iN[nt]]),'-w',lw=1.5)
    ax.plot((lons[jN[nt],i0[nt]],lons[jN[nt],iN[nt]]),(lats[jN[nt],i0[nt]],lats[jN[nt],iN[nt]]),'-w',lw=1.5)

    #ax.text(iloc[nt],jloc[nt],reg_labels[nt],color=[0,0,0])
##############################################
# -@moeba---------------------------------- #
##############################################

# FILE DETAILS
dhw_file = '/data/external/P1/Data/CORTAD/Version4/cortadv4_TSA_DHW_coral.nc'
fid = nc.Dataset(dhw_file)

cmap_file = '/home/frederic/python/cmap/dhw_noaa.cmap'

shp_file = '/t3/workdir/liz/external_data/CT_ATLAS/Coral_Triangle_Boundary/Coral_Triangle_Boundary_Line.shp'
sfid = shp.Reader(shp_file)

sstfile = '/data/external/P1/Data/CORTAD/Version4/cortadv4_FilledSST_coral.nc'
sstfid = nc.Dataset(sstfile)

# EXTRACT VARIABLES
cor_time = sstfid.variables['time'][:]
ref = dt.datetime(1900,1,1,0,0)
YOI = 1998
yrs = np.arange(1982,YOI+1)
xnew = range(1,365+1)

#itime = range(1470,1521+1) # 2010 
#itime = range(844,895+1) # 1998
itime = index_weeks(YOI)
time = fid.variables['time'][itime]
plot_time = convert_time(time)
lats = fid.variables['lat'][:]
lons = fid.variables['lon'][:]

DHW = fid.variables['dhw'][itime,:,:]
DHW = np.ma.masked_where(DHW<0,DHW)
DHW_max = np.max(DHW,axis=0)

# Edit colormap
cmap = np.loadtxt(cmap_file)
lev_seven = np.array((159,212,0))
cmap2 = np.vstack((cmap[0:7,:],lev_seven,cmap[7:,:]))
cmap = cmap2/256.
dhw_noaa = colors.ListedColormap(cmap)

# dhw limits
cmin=0; cmax=16

# FIGURE DETAILS
fig = plt.figure(figsize=(15,9))
gs1 = gridspec.GridSpec(4, 3)
ax1 = fig.add_subplot(gs1[:])

# main map domain
lats1 = [-22.,25.]
lons1 = [90.,170.]
# inset map domain
lats2 = [12,15]
lons2 = [119.5,122.5]

# REGION BOUNDS
j0 = [700, 770, 940, 755, 670,  860]
jN = [1050, 820, 990, 805, 720,  875]
i0 = [460,530, 640, 790, 970,  735]
iN = [1100, 580, 690, 840, 1020, 750]
#reg_labels = ['all','i','ii','iii','iv','v']
# PRIMARY CT MAP
m = Basemap(llcrnrlon=lons1[0],llcrnrlat=lats1[0],urcrnrlon=lons1[1],urcrnrlat=lats1[1],resolution='i',ax=ax1)
cs = m.pcolormesh(lons,lats,DHW_max,vmin=cmin,vmax=cmax,cmap=dhw_noaa,latlon=True)
polygon_patch(m,ax1)
plt_shp_coord(sfid,m)

m.drawmeridians([100,120,140,160], labels=[0,0,1,0], fmt='%d', fontsize=14)
m.drawparallels([-20,0,20], labels=[1,0,0,0], fmt='%d', fontsize=14)

cb_ticks = np.linspace(cmin,cmax,5)
clb = m.colorbar(cs,location='bottom',ticks=[])
#clb.set_label('Maximum Degree Heating Week Signal', fontsize=14, labelpad=25)

#for t in clb.ax.get_yticklabels():
#    t.set_fontsize(18)

# SST Clim plots
fig2 = plt.figure()
# gs2 = gridspec.GridSpec(len(j0), 1)
#gs2.update(left=pos.x1-.05, right = pos.x1+.12,hspace=0)

mon_vals = np.array(range(3,12,4))
months = pltd.MonthLocator(bymonth=mon_vals, bymonthday=1)#interval=3) 
month_labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
labels = [month_labels[nl-1] for nl in mon_vals]

ymins = [23.5,25, 22, 25.5, 26, 24.5]
ymaxs = [31.5,32, 33, 32.5, 32, 31.5]

yt1 =   [24,26, 24, 27, 27, 26]
yt3 =   [30,31, 31, 31, 31, 30]

iloc = [lons[jN[0],iN[0]],lons[jN[1],iN[1]],lons[jN[0],iN[0]],lons[jN[1],iN[1]],lons[j0[2],i0[2]]-1,lons[j0[3],i0[3]]-1,lons[jN[4],iN[4]]]
jloc = [lats[jN[0],iN[0]],lats[jN[1],iN[1]],lats[jN[0],iN[0]],lats[jN[1],iN[1]],lats[j0[2],i0[2]]-1,lats[j0[3],i0[3]]-1,lats[jN[4],iN[4]]]

#for nt in range(len(j0)):
for nt in range(5,6):
    ax = fig2.add_subplot(111)
    plt_sst_ts(nt,ax)
    if nt == len(j0)-1:
       axins = plt_inset(ax1)
       plt_box(nt,axins)

    else:
       plt_box(nt,ax1)

#gs2.update(top=gs1.top, bottom=gs1.bottom)
#fig.text(0.94, 0.45, r'SST ($^{o}$C)', rotation=-90, va="center",fontsize=14)
plt.show()
