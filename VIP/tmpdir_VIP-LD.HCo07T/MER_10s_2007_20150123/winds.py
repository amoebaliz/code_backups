import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap


# MONTH NAMES
title_val = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

fig = plt.figure(num = 1, figsize = (25,10))
mngr = plt.get_current_fig_manager()
geom = mngr.window.geometry()
x,y,dx,dy = geom.getRect()
mngr.window.setGeometry(200,150,dx,dy)
# QUIVER ARROW FREQUENCY
afreq = 10
# ITERATE OVER AND DISPLAY ALL MONTHS
for J in np.arange(len(title_val)):
#for J in [3]:
     # FILE NAME
     ncfile_wind = title_val[J] + '_wind.nc'
     ncfile_current = title_val[J] + '_current.nc'

     fid = nc.Dataset(ncfile_wind, 'r')
     fid2 = nc.Dataset(ncfile_current, 'r')

     Uwind = fid.variables['Uwind'[:]]
     Vwind = fid.variables['Vwind'[:]]
     zeta = fid.variables['zeta'[:]]

     U = fid2.variables['u'[:]]	     
     V = fid2.variables['v'[:]]
     u = np.mean(U[:,-1,:,:],axis=0)

     Uwind_avg = np.mean(Uwind[:,:,:], axis = 0)
     Vwind_avg = np.mean(Vwind[:,:,:], axis = 0)
     Zeta_avg = np.mean(zeta[:,1:-1,1:-1], axis = 0)
     # INTERPOLATE TO RHO-POINTS 
     U = np.squeeze((U[:,-1,1:-1,0:-1]+U[:,-1,1:-1,1:])/2)
     V = np.squeeze((V[:,-1,0:-1,1:-1]+V[:,-1,1:,1:-1])/2)

     U_avg = np.mean(U, axis = 0)
     V_avg = np.mean(V, axis = 0)
     # CALCULATE CURRENT SPEED 
     current_speed = np.sqrt(np.square(U_avg)+np.square(V_avg))	

     fid.close()
     fid2.close()

     x = np.arange(0,U_avg.shape[1])
     y = np.arange(0,V_avg.shape[0])

     X, Y = np.meshgrid(x, y)

     if J > 0:
          Q1.remove()
          Q2.remove()
     plt.subplot(131)
     Q1 = plt.quiver(Uwind_avg[::afreq,::afreq],Vwind_avg[::afreq,::afreq])
     plt.title(title_val[J] + ' winds')     


     plt.subplot(132)
     #C1 = plt.contourf(X,Y,current_speed,100)
     C1 = plt.contourf(X,Y,U_avg,100)
     #plt.colorbar()
     Q2 = plt.quiver(X[::afreq,::afreq],Y[::afreq,::afreq],U_avg[::afreq,::afreq],V_avg[::afreq,::afreq])
     plt.title(title_val[J] + ' currents')
     plt.tight_layout(pad=4, w_pad=5, h_pad=1.0)
 
     plt.subplot(133)
     C2 = plt.contourf(X,Y,Zeta_avg,100)
     plt.title(title_val[J] + ' SSH')
     #plt.colorbar()
     #plt.clim([.5,1])
     plt.pause(1)

plt.show()


# MAKING MONTHLY FILES:
# ncrcat -v Uwind,Vwind,zeta  -n 31,4,1 VIP-LD.HCo07T_his_0001.nc JAN_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 28,4,1 VIP-LD.HCo07T_his_0032.nc FEB_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 31,4,1 VIP-LD.HCo07T_his_0060.nc MAR_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 30,4,1 VIP-LD.HCo07T_his_0091.nc APR_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 31,4,1 VIP-LD.HCo07T_his_0121.nc MAY_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 30,4,1 VIP-LD.HCo07T_his_0152.nc JUN_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 31,4,1 VIP-LD.HCo07T_his_0182.nc JUL_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 31,4,1 VIP-LD.HCo07T_his_0212.nc AUG_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 30,4,1 VIP-LD.HCo07T_his_0243.nc SEP_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 31,4,1 VIP-LD.HCo07T_his_0273.nc OCT_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 30,4,1 VIP-LD.HCo07T_his_0304.nc NOV_wind.nc
# ncrcat -v Uwind,Vwind,zeta  -n 28,4,1 VIP-LD.HCo07T_his_0335.nc DEC_wind.nc


# ncrcat -v zeta,u,v -n 31,4,1 VIP-LD.HCo07T_avg_0001.nc JAN_current.nc
# ncrcat -v zeta,u,v -n 28,4,1 VIP-LD.HCo07T_avg_0032.nc FEB_current.nc
# ncrcat -v zeta,u,v -n 31,4,1 VIP-LD.HCo07T_avg_0060.nc MAR_current.nc
# ncrcat -v zeta,u,v -n 30,4,1 VIP-LD.HCo07T_avg_0091.nc APR_current.nc
# ncrcat -v zeta,u,v -n 31,4,1 VIP-LD.HCo07T_avg_0121.nc MAY_current.nc
# ncrcat -v zeta,u,v -n 30,4,1 VIP-LD.HCo07T_avg_0152.nc JUN_current.nc
# ncrcat -v zeta,u,v -n 31,4,1 VIP-LD.HCo07T_avg_0182.nc JUL_current.nc
# ncrcat -v zeta,u,v -n 31,4,1 VIP-LD.HCo07T_avg_0212.nc AUG_current.nc
# ncrcat -v zeta,u,v -n 30,4,1 VIP-LD.HCo07T_avg_0243.nc SEP_current.nc
# ncrcat -v zeta,u,v -n 31,4,1 VIP-LD.HCo07T_avg_0273.nc OCT_current.nc
# ncrcat -v zeta,u,v -n 30,4,1 VIP-LD.HCo07T_avg_0304.nc NOV_current.nc
# ncrcat -v zeta,u,v -n 27,4,1 VIP-LD.HCo07T_avg_0335.nc DEC_current.nc

# ncrcat -v temp,salt,rho -n 31,4,1 VIP-LD.HCo07T_avg_0001.nc JAN_chem.nc
# ncrcat -v temp,salt,rho -n 28,4,1 VIP-LD.HCo07T_avg_0032.nc FEB_chem.nc
# ncrcat -v temp,salt,rho -n 31,4,1 VIP-LD.HCo07T_avg_0060.nc MAR_chem.nc
# ncrcat -v temp,salt,rho -n 30,4,1 VIP-LD.HCo07T_avg_0091.nc APR_chem.nc
# ncrcat -v temp,salt,rho -n 31,4,1 VIP-LD.HCo07T_avg_0121.nc MAY_chem.nc
# ncrcat -v temp,salt,rho -n 30,4,1 VIP-LD.HCo07T_avg_0152.nc JUN_chem.nc
# ncrcat -v temp,salt,rho -n 31,4,1 VIP-LD.HCo07T_avg_0182.nc JUL_chem.nc
# ncrcat -v temp,salt,rho -n 31,4,1 VIP-LD.HCo07T_avg_0212.nc AUG_chem.nc
# ncrcat -v temp,salt,rho -n 30,4,1 VIP-LD.HCo07T_avg_0243.nc SEP_chem.nc
# ncrcat -v temp,salt,rho -n 31,4,1 VIP-LD.HCo07T_avg_0273.nc OCT_chem.nc
# ncrcat -v temp,salt,rho -n 30,4,1 VIP-LD.HCo07T_avg_0304.nc NOV_chem.nc
# ncrcat -v temp,salt,rho -n 28,4,1 VIP-LD.HCo07T_avg_0335.nc DEC_chem.nc
