import netCDF4 as nc
import numpy as np
import pyroms
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime as dt

def get_data(file,bounds):
    fid = nc.Dataset(file)
    ssh = np.squeeze(np.mean(fid.variables['zeta'][bounds[0]:bounds[1]+1,:,:],axis=0))
    u = np.squeeze(np.mean(fid.variables['u'][bounds[0]:bounds[1]+1,:,:,:],axis=0))
    u_avg = (u[:,1:] + u[:,:-1])/2.
    v = np.squeeze(np.mean(fid.variables['v'][bounds[0]:bounds[1]+1,:,:,:],axis=0))
    v_avg = (v[1:,:] + v[:-1,:])/2.

    return ssh, u_avg, v_avg

def plot_map_vals(ssh,u,v):
    fig = plt.figure(figsize=(9,6))
    m = Basemap(llcrnrlon=117,llcrnrlat=5,urcrnrlon=126,urcrnrlat=14,resolution='h')
    m.drawcoastlines()
    m.fillcontinents()
    pcol = m.pcolor(lons, lats, ssh, vmin=.40,vmax=1.1) 
    cbar = plt.colorbar(pcol)
    Q = m.quiver(lons[::afreq,::afreq],lats[::afreq,::afreq],u[::afreq,::afreq],v[::afreq,::afreq],scale = 1/.3, width = 0.003)

CORAL = pyroms.grid.get_ROMS_grid('CORAL')

ct_Ilats = [310,515]
ct_Ilons = [311,564]
lats = CORAL.hgrid.lat_rho[ct_Ilats[0]:ct_Ilats[1]+1,ct_Ilons[0]:ct_Ilons[1]+1]
lons = CORAL.hgrid.lon_rho[ct_Ilats[0]:ct_Ilats[1]+1,ct_Ilons[0]:ct_Ilons[1]+1]
# FILE DETAILS
fil_dir = '/data/external/P4/workdir/liz/ncfile_kitchen/ctroms_vip/'
#time_vals = range(60,243)
afreq = 5
yrs = (2005,2006,2006,2006)
bounds = np.array([[304,333],[0,30],[151,180],[212,242]])
n=0
for nyr in yrs:

    file = fil_dir + 'ct_vip_zetavel_' + str(nyr) + '_nointerp.nc'
    ssh, u, v = get_data(file,bounds[n,:])
    plot_map_vals(ssh,u,v)
    n+=1 
plt.show()    
