import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt

# DEFINITIONS
def plot_square(ax,pt1,pt2):
    x_vals = [pt1[0],pt1[0],pt2[0],pt2[0],pt1[0]] 
    y_vals = [pt1[1],pt2[1],pt2[1],pt1[1],pt1[1]]
    for nt in range(5):
        ax.plot(x_vals[nt:(nt+2)],y_vals[nt:(nt+2)],'-k')

def subset_bounds(pt1,pt2):
    x1 = np.min([pt1[0],pt2[0]])
    x2 = np.max([pt1[0],pt2[0]])
    y1 = np.min([pt1[1],pt2[1]])
    y2 = np.max([pt1[1],pt2[1]])
    return x1, x2, y1, y2

def convert_time(time):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(days=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

def plot_map(bathy_file):
    fig,ax = plt.subplots(1)
    fid = nc.Dataset(bathy_file)
    bathy = fid.variables['h'][:]
    mask_vals = fid.variables['mask_rho'][:]
    masked_bathy = np.ma.masked_where(mask_vals==0,bathy)
    C = ax.pcolor(masked_bathy)
    plt.colorbar(C)
    ax.set_xlim(0,masked_bathy.shape[1])
    ax.set_ylim(0,masked_bathy.shape[0])
    return ax

def plot_model_ts(pt1,pt2):
    fig,ax = plt.subplots(1, figsize = (16,5))
    x1,x2,y1,y2 = subset_bounds(pt1,pt2)
    for nmod in range(len(models)):
        fid = nc.Dataset(models[nmod])
        time = fid.variables[times[nmod]][:]
        plot_times = convert_time(time)
        temp = np.squeeze(fid.variables[temps[nmod]][:])
        plot_data = np.mean(np.mean(temp[:,y1:y2,x1:x2],axis=2),axis=1)
        ax.plot(plot_data,color_vals[nmod])
#        ax.plot(plot_times,plot_data,color_vals[nmod]) 
#    ax.xaxis_date()

# SET UP
bnds1 = np.array([[150,107],[369,208],[314,164]])
bnds2 = np.array([[200,138],[408,227],[329,179]])

vip_trans_file = '/t3/workdir/liz/scripts/VIP_analyses/vip_u_transport_col278_rows144-178.nc'

#models = ['CORTAD_VIP_sst_1996-1999.nc',\
#          'CTROMS_VIP_sst_1996-1999.nc',\
#          'VIP_sst_1996-1999.nc']

models = ['CTROMS_VIP_sst_1996-1999.nc',\
          'VIP_sst_1996-1999.nc']

times = ['ocean_time','ocean_time']

temps = ['temp','temp']

color_vals = ['-r','-b']

label_vals = ['CORTAD','CT_ROMS','VIP']

vip_bathy_file = '/t3/workdir/liz/VIP/Inputs/Grid/VIP_grd_high_res_bathy_interp.nc'
ax2 = plot_map(vip_bathy_file)

for ndom in range(bnds1.shape[0]):
    pt1 = bnds1[ndom,:]
    pt2 = bnds2[ndom,:]
    plot_model_ts(pt1,pt2)
    plot_square(ax2,pt1,pt2)

plt.show()
