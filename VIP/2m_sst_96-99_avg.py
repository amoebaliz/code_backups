import numpy as np
import netCDF4 as nc
import pyroms
import matplotlib.pyplot as plt
import time
import datetime as dt

def get_zw(grd,zeta):
    h = grd.vgrid.h
    hc = grd.vgrid.hc
    Np = grd.vgrid.Np
    s_w = grd.vgrid.s_w
    Cs_w = grd.vgrid.Cs_w
    Vtrans = grd.vgrid.Vtrans
    z = pyroms.vgrid.z_r(h, hc, Np, s_w, Cs_w, zeta, Vtrans) 
    z = z[0:]
    return z

# ADAPTED FROM PYROMS ZLAYER
def zlayer(var, zeta, grd, h2, Cpos='rho', vert=False):
    # compute the depth on Arakawa-C grid position
    # ASSUMED: 
    # 1) Cpos = rho
    # 2) vert == False (default)
    # 3) NOT a Spherical grid
    
    # hand the ZETA field to the z_w function for given time point
    z = get_zw(grd,zeta)
    mask = grd.hgrid.mask_rho[:]

    #set var to zero where tra is masked for the sum (??)
    mask = np.tile(mask, (var.shape[0],1,1))
    var = np.where(mask == 1, var, 0)
    zlayer = np.zeros((var.shape[1], var.shape[2]))
    # layer thicknesses
    dz = z[1:,:,:] - z[:-1,:,:]
    # iterating all lat lons
    for jj in range(var.shape[1]):
        for ji in range(var.shape[2]):
            # initialize ratio array; len = #of depth levels
            ratio = np.ones(var.shape[0])
            # CRITERIA: uppermost depth level must be above bottom depth of bottom layer
            # if h1 is not None:
            #   if np.abs(h1) < np.abs(z[0,jj,ji]):
	    #      idx1 = np.where(z[:,jj,ji] > -np.abs(h1))
            #      if np.any(idx1):
                     # which depth is on the cusp of the layer (i.e., only a fraction is part of the layer)
            #         idx1 = idx1[0][0]
            #         r = (-np.abs(h1) - z[idx1-1,jj,ji]) / \
            #             (z[idx1,jj,ji] -z[idx1-1,jj,ji])
            #         ratio[idx1:] = 0.
            #         ratio[idx1-1] = r
            #   else:
            #         ratio[:] = 0.
             
            if h2 is not None:
               sum_dz = np.cumsum(dz[::-1,jj,ji])                     
               sum_dz = sum_dz[::-1]
               idx2 = np.where(sum_dz > np.abs(h2))
                
               if np.any(idx2):
                  idx2 = idx2[0][-1]

                  # NOT top layer 
                  if idx2 < (len(sum_dz)-1):
                     r = np.float((abs(h2) - sum_dz[idx2+1]) / \
                         (sum_dz[idx2] - sum_dz[idx2+1]))
                  # TOP LAYER
                  else:
                     r = np.float(abs(h2) / \
                         sum_dz[-1]) 
                  ratio[:idx2] = 0
                  ratio[idx2] = r

            # WEIGHTED AVERAGE
               zlayer[jj,ji] = np.sum(var[:,jj,ji] * dz[:,jj,ji] * ratio[:])/ \
                               np.sum(dz[:,jj,ji] * ratio[:])    

    return zlayer



#############   CODE   #################
grd = pyroms.grid.get_ROMS_grid('VIP')
h1=None      # Surface 
h2=-2        # Depth (m)

n=0
file_val = 13882 
days = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
data_dir = '/t3/workdir/liz/VIP/Runs/VIP-LD.HCo10T/outputs/'
#data_dir = '/data/external/P4/workdir/liz/ncfile_kitchen/VIP_Interpolations/ct_temp/'
#zeta_dir = '/t3/workdir/liz/scripts/VIP_analyses/CTROMS_VIP_zeta_1996-1999.nc'
for nyr in range(1996,1996+1):
    if nyr%4 == 0:
       days[1] = 29
    else: days[1] = 28

    #for nmon in range(12):
    for nmon in [0]:
        #for nday in range(days[nmon]):
        for nday in [0]:
            file = data_dir + str(nyr)+'/VIP-LD.HCo10T_avg_' + str(nyr) + '-' + \
            str(nmon+1).zfill(2) + '-' + str(nday+1).zfill(2) + 'T00:00:00.nc'
            #temp_file = data_dir + 'coral_avg_' + str(file_val) + '_VIP.nc'
            #zeta_file = zeta_dir
            fid = nc.Dataset(file)
            #fid2 = nc.Dataset(zeta_file)
            temp = np.squeeze(fid.variables['temp'][:])
            print temp.shape
            zeta = np.squeeze(fid.variables['zeta'][:])
            tx = fid.variables['ocean_time']
            avg2m_temp = zlayer(temp,zeta,grd,h2)
           
            fid_test = nc.Dataset('/t3/workdir/liz/VIP/Runs/VIP-LD.HCo10T/outputs/1996/VIP-LD.HCo10T_avg_1996-01-01T00:00:00.nc')
            sst_test = np.squeeze(fid_test.variables['temp'][:,-1,:,:])
 
            plot_sst = np.ma.masked_where(np.isnan(avg2m_temp),avg2m_temp) 
            plt.figure()
            plt.pcolor(plot_sst-sst_test,vmin=0, vmax=0.03)
            plt.colorbar()
            plt.show()
 
            # Write time and sst to netcdf file 
            newfile = 'sst_2m_temp_vip' + str(n).zfill(4) + '.nc'
            nfid = nc.Dataset(newfile, 'w')
            nfid.createDimension('ocean_time', None)
            nfid.createDimension('eta_rho', temp.shape[-2])
            nfid.createDimension('xi_rho', temp.shape[-1])
            time = nfid.createVariable('ocean_time','f8',('ocean_time',))
            time[:] = tx[:]
            time.units = tx.units

            ref = ref = dt.datetime(1900,1,1,0,0)
            newdate = ref + dt.timedelta(seconds=(time[0]))
            print newdate

            sst = nfid.createVariable('layer_temp','f4',('ocean_time','eta_rho','xi_rho',))
            sst[0,:,:] = avg2m_temp
            sst.long_name = 'average temp in top 2m' 
            nfid.close()
            
            print n
            n+=1
            file_val+=1
