import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import pyroms as py

# MONTH NAMES
title_val = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

VIP = py.grid.get_ROMS_grid('VIP')

nc_grid_file = '/t1/scratch/liz/Inputs/VIP-LD.HCo11T/Grid/VIP_grd_high_res_bathy_interp2.nc'
fid = nc.Dataset(nc_grid_file, 'r')
mask_vals = fid.variables['mask_rho'[:]]
h = fid.variables['h'[:]]
s_rho = fid.variables['s_rho'[:]]

mask_grid = mask_vals[:]
depth = h[1:-1,1:-1]
s_rho = s_rho[:]

fid.close()

class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
        if len(self.xs)==2:
           self.line.figure.canvas.mpl_disconnect(self.cid)


#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.set_title('click to build transect')
#line, = ax.plot([], [])  # empty line
#linebuilder = LineBuilder(line)
#plt.pcolor(mask_grid)

#plt.draw()
#plt.pause(6)
#plt.close()
# CALCULATE stuff

#min_row = 150
#max_row = 205
#rows = np.arange(min_row,max_row+1)
#col = 320

#min_row = 140
#max_row = 180
#rows = np.arange(min_row,max_row+1)
#col = 280

#min_row = 50
#max_row = 145
#rows = np.arange(min_row,max_row+1)
#col = 188

row = 162
max_col = 300
min_col = 240
cols = np.arange(min_col,max_col+1)

#local_seg_width = (depth[rows,col]/len(s_rho)).reshape(1,len(rows))
#dep_vec = (-1*(np.arange(len(s_rho))+.5)).reshape((len(s_rho),1))
#local_dep_vals = np.multiply(local_seg_width,dep_vec)
X = np.tile(cols,(len(s_rho),1))
# TRANSECT LOCATION PLOT
fig1 = plt.figure()
plt.pcolor(mask_grid)
plt.plot([min_col,max_col],[row,row],'g^-')
plt.tight_layout(pad=4, w_pad=5, h_pad=1.0)
plt.xlim([0,mask_grid.shape[1]])
plt.ylim([0,mask_grid.shape[0]])

# FIGURE CONFIGURATION

fig2 = plt.figure(num = 2, figsize = (25,7))
mngr = plt.get_current_fig_manager()
geom = mngr.window.geometry()
x,y,dx,dy = geom.getRect()
mngr.window.setGeometry(200,50,dx,dy)

n = 0
for J in np.arange(len(title_val)):
#for J in [3]:
     # FILE NAME
     ncfile_current = title_val[J] + '_current.nc'
     #ncfile_chem = title_val[J] + '_chem.nc'
     fid = nc.Dataset(ncfile_current,'r')
     #fid2 = nc.Dataset(ncfile_chem,'r')

     #u = fid.variables['u'[:]]
     #v = fid.variables['v'[:]]
     temp = fid2.variables['temp'[:]]
     #u = u[:]
     temp = temp[:]
     #u_avg = np.squeeze(np.mean(u[:,:,:,:],axis=0))
     temp_avg = np.squeeze(np.mean(temp[:,:,:,:],axis=0))
     #u_slice,u_Z,u_X,u_Y = py.tools.islice(u_avg,col,VIP,'u')
     temp_slice, temp_Z, temp_X, temp_Y = py.tools.jslice(temp_avg,row,VIP,'rho')

     n = n+1
     plt.subplot(3,4,n)
     C1 = plt.contourf(X,temp_Z[:,cols],temp_slice[:,cols],100) 
  
     plt.clim(13,27)
     plt.title(title_val[J] + ' Temperature (oC)')
     if n > 8:     
         plt.xlabel('J Direction')

plt.show()


