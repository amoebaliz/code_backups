# calculating transport index
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms.tools as pyt
import pyroms as py
import datetime as dt
import matplotlib.dates as pltd
import matplotlib.animation as animation
import matplotlib.gridspec as gsp

vip_dir = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/'
#ctroms_interp_dir = ''

ref = dt.datetime(1900,1,1,0,0)

# TRANSECT POINTS [i,j]
vip1 = [278,144]; vip2 = [278,178]

# GRIDS
VIP = py.grid.get_ROMS_grid('VIP')

afreq = 10
col = [190,242,280,320]
min_row = [70, 100, 140, 140]
max_row = [150, 150, 180, 200]
let_val = ['A','B','C','D']



vip_nums = np.arange(5,7 + 1)
vip_nums = np.arange(5, 222 + 1)
transect_store = np.ma.zeros((len(vip_nums),50,35))
vip_transport = np.zeros(len(vip_nums))
time = np.zeros(len(vip_nums))
n = 0
for nt in vip_nums:
   vip_file = vip_dir + 'VIP-LD.HCo07T_avg_' + str(nt).zfill(4) + '.nc'
   vip_fid = nc.Dataset(vip_file,'r')
   u = vip_fid.variables['u'][:]; u = np.squeeze(u)
   # temp = vip_fid.variables['temp'][:]; temp = np.squeeze(temp)
   time[nt-vip_nums[0]] = vip_fid.variables['ocean_time'][:]
   s_rho = vip_fid.variables['s_rho'][:];vs_rho = s_rho[:]
   transect_store[nt-vip_nums[0],:,:], z, lon, lat = pyt.transect(u,vip1[0],vip2[0],vip1[1],vip2[1],VIP,'u',spval=1e+37)
      
   
fid2 = nc.Dataset('vip_temp_transect.nc','w')
fid2.createDimension('ocean_time', None)
fid2.createDimension('depth', 50)
fid2.createDimension('lateral', vip2[1]-vip1[1]+1)

X = np.tile(np.arange(0,z.shape[1]),(z.shape[0],1))
X_var = fid2.createVariable('X','f8',('depth','lateral',))
X_var[:] = X
ocean_time = fid2.createVariable('ocean_time','f8',('ocean_time',))
ocean_time.long_name = 'averaged time since initialization'
ocean_time.units = 'seconds since 1900-01-01 00:00:00'
ocean_time[:] = time

z_field = fid2.createVariable('z_field','f8',('depth','lateral',))
z_field.units = 'meter'
z_field[:] = z

temp = fid2.createVariable('temp','f8',('ocean_time','depth','lateral',))
temp.long_name = 'time-averaged potential temperature' 
temp.units = 'Celsius'
temp[:] = transect_store
fid2.close()

######## PLOTTING ANIMATION
# contour levels

clev = 1
levs = np.linspace(-1*clev,clev,100,endpoint=True)
X = np.tile(np.arange(0,z.shape[1]),(len(s_rho),1))

fig,axs = plt.subplots(figsize=(4,4))
#gs = gsp.GridSpec(1,4,width_ratios = np.array(max_row)-np.array(min_row)+1)
#gs.update(bottom = 0.2, left = 0.1, right = 0.88, top = 0.87, wspace=0.2)

axs.set_axis_bgcolor((0.7,0.7,0.7))
cbar_ax = fig.add_axes([0.93,0.22,0.015,0.5])
cmap = plt.cm.bwr
cmap.set_bad('k',1.)

def animate(i):
    cont = axs.contourf(X,z,np.squeeze(transect_store[i,:,:]),levs,cmap=cmap)
    axs.set_title('day %i' % i)
    fig.colorbar(cont,cax=cbar_ax,ticks = [-1*(clev),0,clev])
    return cont

ani = animation.FuncAnimation(fig, animate, frames=len(vip_nums))
#mywriter = animation.FFMpegWriter()
#ani.save('animation.mp4', writer=mywriter)
plt.show()
