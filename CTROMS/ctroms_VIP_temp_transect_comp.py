# COMPARING TEMPERATURE OVER DEPTHS
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import pyroms as py
import datetime as dt

# INDICIES INFORMATION
# note: the first 4 average files have are empty so nplus must be >4
nplus = [10, 15, 20, 25, 30, 35]
col = [100, 190, 240, 280, 320, 400] 
min_row = [87, 70, 100, 145, 145, 150] ; max_row = [335, 150, 153, 180, 200, 335]
tran_id = ['A','B','C','D','E','F']
ref = dt.datetime(1900,1,1,0,0)
# GRID INFORMATION
VIP = py.grid.get_ROMS_grid('VIP')
mask_vals = VIP.hgrid.mask_rho
s_rho = VIP.vgrid.s_rho

mask_grid = mask_vals[:]
s_rho = s_rho[:]

vip_init = 1
ctr_init = 17900
levels = np.linspace(0,35,num=100,endpoint=True)
# MAP FOR TRANSECT LOCATIONS
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
plt.pcolor(mask_grid)
plt.xlim([0,mask_grid.shape[1]])
plt.ylim([0,mask_grid.shape[0]])
ax1.set_xlabel('J Direction',fontsize=14,labelpad=10)
ax1.set_ylabel('I Direction',fontsize=14,labelpad=12)
for nt in range(len(col)):
# for nt in range(2,3): 
   print nt  
   rows = np.arange(min_row[nt],max_row[nt]+1)
   X = np.tile(rows,(len(s_rho),1))
 
   # DATA FILES
   vip_file = '/t3/workdir/liz/MODELS/VIP/PostProc/1996_VIP_average.nc'


#   ctroms_file = '/t3/workdir/liz/MODELS/VIP/PostProc/1996_CT_average.nc'

   fid_vip = nc.Dataset(vip_file,'r')
#   fid_ctroms = nc.Dataset(ctroms_file,'r')

   # DATA VARIABLES
   time_vip = fid_vip.variables['ocean_time'][:]; vtm = time_vip[:]
   temp_vip = fid_vip.variables['temp'][:]; vtemp = temp_vip[:]
   vtemp = np.squeeze(vtemp)
   
#   time_ctroms = fid_ctroms.variables['ocean_time'][:]; ctm = time_ctroms[:]
#   temp_ctroms = fid_ctroms.variables['temp'][:]; ctemp = temp_ctroms[:]
#   ctemp = np.squeeze(ctemp)

   # I SLICE   
   vip_slice, vip_Z, vip_X, vip_Y = py.tools.islice(vtemp,col[nt],VIP,'rho')
   #ctr_slice, ctr_Z, ctr_X, ctr_Y = py.tools.islice(ctemp,col[nt],VIP,'rho')                                                               
   # FIGURE AND PLOTTING
   f, (ax2,ax4) = plt.subplots(1,3,sharey = True,figsize = (16,7))
   f.tight_layout(pad=4,w_pad=None)
   # ORIGINAL VIP MODEL
   a = ax2.contourf(X,vip_Z[:,rows],vip_slice[:,rows],levels)
   a2 = ax2.contour(a, levels=levels[::3],colors='k',hold='on')
   ax2.set_title('VIP-old')
   ax2.set_ylabel('Depth (m)',fontsize=14,labelpad=12)

   # CT-ROMS MODEL 
 #  c = ax4.contourf(X,ctr_Z[:,rows],ctr_slice[:,rows],levels)
  # c2 = ax4.contour(c, levels=levels[::3],colors='k',hold='on')
  # ax4.set_title('CT-ROMS')
  # cbar2 = f.colorbar(c,cmap='jet')
   #plt.colorbar(a, use_gridspec=True)
  # newdate = ref + dt.timedelta(seconds=ctm[0])
   
   title = 'Temperature Transect ' + tran_id[nt] + ' for ' + str(newdate.date())
   plt.suptitle(title,fontsize=14)
   # PRINT TRANSECT LOCATION ON GRID
   ax1.plot([col[nt],col[nt]],[min_row[nt],max_row[nt]],'g^-')
   x = col[nt]+5; y = (min_row[nt]+max_row[nt])/2; s = tran_id[nt]
   ax1.text(x,y,s,color='w') 

plt.show()
