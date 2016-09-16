import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

def get_sst(j):
    fid = nc.Dataset(ncfiles[j])
    sst = np.squeeze(fid.variables[varname[j/2]][:])
    if np.sum(np.isnan(sst))>0:
       sst = np.ma.masked_where(np.isnan(sst),sst)
    return sst

ncfiles = ['sst_-1rho_VIP_96-99_avg.nc', 'sst_-1rho_ctroms_VIP_96-99_avg.nc',\
           'sst_2m_VIP_96-99_avg.nc',    'sst_2m_ctroms_VIP_96-99_avg.nc']

varname = ['temp','layer_temp']

# Create figure
fig,ax = plt.subplots(2,2,sharey = 'row',sharex = 'col')
plt.tight_layout(pad=2.0, w_pad=0.5, h_pad=2.0)
fig.subplots_adjust(right=0.85,top=0.93)
cbar_ax = fig.add_axes([0.88,0.2,0.02,.65]) 

# Draw Pcolor Plots
for nt in range(4):
#for nt in range(len(ncfiles)):
    a = nt/2 
    b = nt-2*(nt/2)
    sst = get_sst(nt)
    C = ax[a,b].pcolor(sst,vmin=27,vmax=30)
    ax[a,b].set_xlim(0,sst.shape[1])
    ax[a,b].set_ylim(0,sst.shape[0])
    ax[a,b].set_axis_bgcolor((0.5, 0.5, 0.5))
plt.colorbar(C,cax=cbar_ax)
plt.show()
