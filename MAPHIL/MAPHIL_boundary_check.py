import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import pyroms

def convert_time(time,i):
    # NOTE: TIME VARIABLES MUST BE IN DAYS SINCES (1900,1,1,0,0)
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        if i > 1:
           day_time = ref + dt.timedelta(seconds=np.float(time[nt]))
        else:
           day_time = ref + dt.timedelta(days=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals 

# ------------------------------------------------------ #
GRD = pyroms.grid.get_ROMS_grid('MaPhil')
#GRD = pyroms.grid.get_ROMS_grid('VIP')
s_rho = 0
eta_rho = 112
#eta_rho = 90
xi_rho = 0
z = GRD.vgrid.z_r
depths = z[0][:,eta_rho,xi_rho]

ncfile = ['/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Boundary/hycom_BDRY.nc',\
          '/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Clim/hycom_clim.nc', \
#       '/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Boundary/hycom_BDRY.nc',\
#        '/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Clim/hycom_clim.nc',\
#        '/t3/workdir/liz/MODELS/MAPHIL/PostProc/OBCFAC_edit/surf_his_temp.nc',\
#        '/t3/workdir/liz/MODELS/MAPHIL/Runs/1_new_inputs/2009/surf_temp.nc',\
#        '/t3/workdir/liz/MODELS/MAPHIL/Runs/1_new_inputs/2009/bot_temp.nc']
#        '/t3/workdir/liz/MODELS/MAPHIL/PostProc/OBCFAC_edit/bot_his_temp.nc']
          '/t3/workdir/liz/MODELS/MAPHIL/PostProc/surf_his_temp.nc',\
          '/t3/workdir/liz/MODELS/MAPHIL/PostProc/bot_his_temp.nc']
#         '/t1/scratch/liz/2_old_tmpdir_MaPhil-LD.HCo05T/2009/surf_his_temp.nc',\
#         '/t1/scratch/liz/2_old_tmpdir_MaPhil-LD.HCo05T/2009/bot_his_temp.nc']
#        '/t3/workdir/liz/MODELS/MAPHIL/PostProc/surf_his_u.nc']#,\
#        '/t3/workdir/liz/MODELS/MAPHIL/PostProc/bot_his_u.nc']
#ncfile = ['/t1/scratch/liz/Inputs/VIP-LD.HCo10T/Boundary/VIP_BRY_y1995-2000.nc',\
#        '/t1/scratch/liz/Inputs/VIP-LD.HCo10T/Clim/CORAL_VIP_CLIM_1995-1999.nc',\
#        '/t3/workdir/liz/MODELS/VIP/PostProc/bot_temp.nc']


var_name = ['temp_west','temp','temp','temp']
#var_name = ['u_west','u','u','u']

#legend_vals = ['BDRY','CLIM','VIP_BOT']
legend_vals = ['BDRY_BOT','CLIM_BOT','MAPHIL_SURF','MAPHIL_BOT','BDRY_SURF','CLIM_SURF']

fig = plt.figure(figsize = (12,10))
ax = fig.add_subplot(111)

for nt in range(len(ncfile)):
    print nt
    fid = nc.Dataset(ncfile[nt])
    time = fid.variables['ocean_time'][:]
    date_vals = convert_time(time,nt) 
    var = np.squeeze(fid.variables[var_name[nt]][:])
    if nt == 0:
       #temp = np.squeeze(temp[:,s_rho,eta_rho])
       ax.plot(date_vals,np.squeeze(var[:,s_rho,eta_rho]),label=legend_vals[nt])
       ax.plot(date_vals,np.squeeze(var[:,49,eta_rho]),label=legend_vals[-2])
    elif nt == 1:
       #temp = np.squeeze(temp[:,s_rho,eta_rho,xi_rho])
       ax.plot(date_vals,np.squeeze(var[:,s_rho,eta_rho,xi_rho]),label=legend_vals[nt])
       ax.plot(date_vals,np.squeeze(var[:,49,eta_rho,xi_rho]),label=legend_vals[-1])

    elif nt>1:
       #temp = np.squeeze(temp[:,eta_rho,xi_rho])
       ax.plot(date_vals,np.squeeze(var[:,eta_rho,xi_rho]),label=legend_vals[nt])

#fid = nc.Dataset('/t1/scratch/liz/Inputs/MaPhil-LD.HCo06T/Initial/CT_init_1996_2_02_ic_MaPhil.nc')
fid = nc.Dataset('/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Initial/hycom_init_2009_11_02_ic_MaPhil.nc')
#fid = nc.Dataset('/t1/scratch/liz/Inputs/VIP-LD.HCo10T/Initial/VIP_Ini_fromCORAL_y1995m09md02.nc')
time = fid.variables['ocean_time'][:]
date_vals = convert_time(time,0)
#temp = np.squeeze(fid.variables[var_name[nt]][:,0,eta_rho,xi_rho])
ax.plot(date_vals,np.squeeze(fid.variables[var_name[nt]][:,0,eta_rho,xi_rho]),label = 'INIT_BOT',lw = 0, c = 'y',marker = 'o', ms = 8)
ax.plot(date_vals,np.squeeze(fid.variables[var_name[nt]][:,49,eta_rho,xi_rho]),label = 'INIT_SURF',lw = 0, c = 'y',marker = 'o', ms = 8)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, numpoints=1,loc=10)
ax.xaxis_date()
#ax.set_xlim(pltd.date2num(dt.datetime(1995,9,1)),pltd.date2num(dt.datetime(2000,1,1)))
ax.set_xlim(pltd.date2num(dt.datetime(2009,11,1)),pltd.date2num(dt.datetime(2010,10,31)))
ax.set_ylim(13,32)
#plt.title('Bottom Temp (oC) at western boundary (eta = 112)')
plt.show()

