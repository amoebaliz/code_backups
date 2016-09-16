# LEARDING MLD CODE
import netCDF4 as nc
import numpy as np
import pyroms_toolbox as pyt
import pyroms as py
import matplotlib.pyplot as plt

VIP = py.grid.get_ROMS_grid('VIP')
CORAL = py.grid.get_ROMS_grid('CORAL') 

vip_dir = '/t1/scratch/liz/tmpdir_VIP-LD.HCo07T/MER_10s_2007_20150323/'
ct_dir = '/data/external/P2/ROMS/CORAL/RUN16/'

vip_nums = range(5,361+1)
ct_nums = range(17005,17361+1)
idx_vals = [50,100,200]
cmin = -15 
cmax = 0

def mldd(i,MODEL):
    
    if MODEL == VIP:
         file = vip_dir + 'VIP-LD.HCo07T_avg_' + str(vip_nums[i]).zfill(4) + '.nc'
         fid = nc.Dataset(file,'r')
         dens = np.squeeze(fid.variables['rho'][:]) 
    elif MODEL == CORAL:
         file = ct_dir + 'coral_avg_' + str(ct_nums[i]).zfill(5) + '.nc'
         fid = nc.Dataset(file,'r')
         salt = np.squeeze(fid.variables['salt'][:]) 
         temp = np.squeeze(fid.variables['temp'][:])
         dens = pyt.seawater.dens(salt,temp)  

    mld = mld_adjust_dens(dens,MODEL)
    return mld

def mld_adjust_dens(dens,grd):
  # ADDING ZSLICE TO RETURN SURF DENSITY FROM A PARTICULAR LOCATION
  z = grd.vgrid.z_r[0]
  surf_dens = np.asarray(py.tools.zslice(dens,-5,grd))
  surf_dens = np.squeeze(surf_dens[0,:,:])
  mld_dens = surf_dens + 0.03 #threshold of 0.03 kg/m^3, following de Boyer Montegut et al., 2004
  mld, lon, lat = py.tools.isoslice(z,dens,mld_dens,grd)
  return mld





for nt in idx_vals:
    v_mld = mldd(nt,VIP)
    c_mld = mldd(nt,CORAL) 

    f, axarr = plt.subplots(2,1, figsize=(14,16))
    f.tight_layout(pad=4,w_pad=None,h_pad=None)
    f.subplots_adjust(right=0.85,top=0.93)
    cbar_ax = f.add_axes([0.9,0.15,0.02,0.7])
    title_val = 'Day ' + str(vip_nums[nt])
    a = axarr[0].pcolor(v_mld,vmin=cmin,vmax=cmax)#; axarr[0,0].set_title('CT-ROMS');axarr[0,0].set_ylabel('Western Boundary')
    axarr[0].set_xlim(0,v_mld.shape[1])
    axarr[0].set_ylim(0,v_mld.shape[0])
    axarr[0].set_axis_bgcolor((0.6,0.6,0.6))
    axarr[0].set_title(title_val)

    axarr[1].pcolor(np.squeeze(c_mld[450:495,365:435]),vmin=cmin, vmax=cmax)#; axarr[0,1].set_title('VIP')
    axarr[1].set_axis_bgcolor((0.6,0.6,0.6))
    plt.colorbar(a,cax=cbar_ax)
plt.show()
