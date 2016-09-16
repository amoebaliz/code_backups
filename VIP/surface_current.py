import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gsp
from mpl_toolkits.basemap import Basemap
import pyroms as py

fig, ax = plt.subplots(1, figsize = (20,8))
fig.tight_layout(pad=2)
plt.subplots_adjust(bottom=0.1, right=0.85, top=0.9, left=0.1)
plt.title('Average VIP Surface Currents', y = 1.04, fontsize=16)
cbar_ax = fig.add_axes([0.88,0.12,0.015,0.7])
ax.set_axis_bgcolor((0.7,0.7,0.7))
ax.set_xlabel('J Direction', labelpad=20, fontsize=14)
ax.set_ylabel('I Direction', labelpad=30, fontsize=14)

afreq = 10
col = [190,242,280,320]
min_row = [70, 100, 140, 140]
max_row = [150, 150, 180, 200]
let_val = ['A','B','C','D']
clev = 0.64
levs = np.linspace(0,clev,100,endpoint=True)
#ncfile = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150313/avg_currents_red.nc'
ncfile = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150313/avg_currents_blue.nc'
fid = nc.Dataset(ncfile, 'r')
U = fid.variables['u'][:]
V = fid.variables['v'][:]
# INTERPOLATE TO RHO-POINTS 
u = np.squeeze((U[:,-1,1:-1,0:-1]+U[:,-1,1:-1,1:])/2)
v = np.squeeze((V[:,-1,0:-1,1:-1]+V[:,-1,1:,1:-1])/2)
u2 = U

u = np.float64(u)
v = np.float64(v)
# CALCULATE CURRENT SPEED 
current_speed = np.sqrt(np.square(u)+np.square(v))
fid.close()
print np.max(current_speed)
x = np.arange(0,current_speed.shape[1])
y = np.arange(0,current_speed.shape[0])
X, Y = np.meshgrid(x, y)

C = ax.contourf(X,Y,current_speed,levs)
fig.colorbar(C,cax=cbar_ax,ticks = [0, clev/4, clev/2, 3*clev/4, clev])
cbar_ax.text(-1.2,1.04,'Current Speed \n (m/s)',multialignment='center')

Q = ax.quiver(X[::afreq,::afreq],Y[::afreq,::afreq],u[::afreq,::afreq],v[::afreq,::afreq],width = 0.0015)


for nt in range(0,4):
  ax.plot([col[nt],col[nt]],[min_row[nt],max_row[nt]],'--k',linewidth=6.0)
  ax.plot([col[nt],col[nt]],[min_row[nt],max_row[nt]],'--w',linewidth=3.0)
  ax.text(col[nt]-1,min_row[nt]-10,let_val[nt])

# DEPTH VELOCITY SECTIONS
VIP = py.grid.get_ROMS_grid('VIP')
U_vals = np.squeeze(u2)
clev = 0.71
levs = np.linspace(-1*clev,clev,100,endpoint=True)

fig=plt.figure(figsize=(17,4))
gs = gsp.GridSpec(1,4,width_ratios = np.array(max_row)-np.array(min_row)+1)
gs.update(bottom = 0.2, left = 0.1, right = 0.88, top = 0.87, wspace=0.2)
cbar_ax = fig.add_axes([0.93,0.22,0.015,0.5])
cmap = plt.cm.bwr
cmap.set_bad('k',1.)

for nt in range(0,4):
  rows = np.arange(min_row[nt],max_row[nt]+1)
  X2 = np.tile(rows,(50,1))
  u_slice, u_Z, u_X, u_Y = py.tools.islice(U_vals,col[nt],VIP,'u')
  print np.max(u_slice)
  if nt == 0:
    axs = fig.add_subplot(gs[nt])
    axs.set_ylabel('Depth (m)')
  else: 
    axs = fig.add_subplot(gs[nt],sharey=axs)
    plt.setp(axs.get_yticklabels(),visible=False)

  C = axs.contourf(X2,u_Z[:,rows],u_slice[:,rows],levs,cmap=cmap)
  axs.set_xlim(min_row[nt], max_row[nt])
  axs.set_axis_bgcolor((0.7,0.7,0.7))
  axs.yaxis.tick_left()
  axs.xaxis.tick_bottom()
  axs.xaxis.set_ticks(np.arange(min_row[nt], max_row[nt]+1,10))
  axs.text(min_row[nt]+1.5,-70,let_val[nt])

fig.colorbar(C,cax=cbar_ax,ticks = [-1*(clev),0,clev])
cbar_ax.text(-2.3,1.1,'Current Velocity \n (m/s)',multialignment='center')
fig.suptitle('Depth Transects of Average VIP Zonal Velocity',x=0.5, fontsize=14)
fig.text(.47,0.04,'I Direction',fontsize=14,transform=ax.transAxes)
plt.show()
