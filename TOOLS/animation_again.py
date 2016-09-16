import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import netCDF4 as nc
from numpy import genfromtxt

ilats = np.arange(300,540)
ilons = np.arange(480,800)

months = ['January','February','March','April','May','June',\
'July','August','September','October','November','December']

yrs = np.arange(1960,2007+1)

# FILE DETAILS
file = '/t3/workdir/liz/scripts/CT_roms/ctroms_zeta_monthly.nc'
fid = nc.Dataset(file,'r')
#lat = fid.variables['lat'][ilats]
#lon = fid.variables['lon'][ilons]
#lons,lats = np.meshgrid(lon,lat)
# SOI DATA AND DETAILS
my_data=genfromtxt('data.txt', usecols=np.arange(0,2))
my_data = my_data[110:685+1,:]

soi = my_data[:,1]
soi_date = my_data[:,0]

soi_la = np.copy(soi); soi_la[soi_la>0]= np.nan
soi_el = np.copy(soi); soi_el[soi_el<0]= np.nan

t_dif = 25
width = 1
plot_ind = np.arange(t_dif)

y_max = np.ceil(np.max(np.abs(soi)))
y_min = -1*y_max
# FIGURE DETAILS
fig = plt.figure(figsize=(10,8))
ax1 = plt.subplot2grid((4,1),(0,0),rowspan=3,axisbg = [.2,.2,.2])
ax2 = plt.subplot2grid((4,1),(3,0))

ax1.set_xlim(0,(len(ilons)-1))
ax1.set_ylim(0,(len(ilats)-1))

plt.tight_layout(pad=2.0, w_pad=0.5, h_pad=2.0)
fig.subplots_adjust(right=0.85,top=0.93)
cbar_ax = fig.add_axes([0.88,0.3,0.02,.55])

plt.hold(True)

def get_zeta(i):
    zeta = fid.variables['zeta'][i,ilats,ilons]
    return zeta

#We need to prime the pump, so to speak and create a quadmesh for plt to work with
# ax1.pcolor(get_zeta(0),vmin=0.4,vmax=1.1)
#plt.colorbar(p1,cax=cbar_ax)
#ax2.bar(plot_ind,soi_la[0:25],width,color='b',edgecolor=[0,0,.5])
#ax2.bar(plot_ind,soi_el[0:25],width,color='r',edgecolor=[.5,0,0])

def animate(i):
    #This is where new data is inserted into the plot.
    ax2.cla()
    cbar_ax.cla()    
    print i
    title_str = str(months[i%12]) + ' ' + str(yrs[i/12])

    p1 = ax1.pcolor(get_zeta(i),vmin=0,vmax=1.0)
    ax1.set_title(title_str,horizontalalignment='right')
    plt.colorbar(p1,cax=cbar_ax)   
    cbar_ax.set_title('SSH\n(m)',x = 1.0, y = 1.05)
 
    ax2.bar(plot_ind,soi_la[i:i+25],width,color='b',edgecolor=[0,0,.5])    
    ax2.bar(plot_ind,soi_el[i:i+25],width,color='r',edgecolor=[.5,0,0])
    ax2.plot([t_dif/2.0,t_dif/2.0],[-4,4],'--k',linewidth=2.0)

    ax2.set_ylim(y_min,y_max) 
    ax2.set_yticks([-4,-2,0,2,4])
    ax2.set_xticks([])
    ax2.set_title('SOI')
    ax2.text(t_dif-4, y_max-1.2,r'$El Ni\~no$', fontsize = 14, color = 'r')
    ax2.text(t_dif-4, y_min+.1,r'$La Ni\~na$', fontsize = 14, color = 'b')

anim = animation.FuncAnimation(fig, animate, frames = range(0,576), blit = False)
anim.save('test.gif',writer='imagemagick',fps=2)

#plt.show()
plt.hold(False)
