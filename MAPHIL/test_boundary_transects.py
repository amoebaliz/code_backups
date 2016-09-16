import numpy as np
import netCDF4 as nc
import pyroms
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def updatefig(i):
    global im1,im2,im3,im4,tx,imon,mon_day
    
    if i>ndays-1:
       mon_day = i+1 - np.sum(ndays[:imon])
    else:
       mon_day = i+1
    if mon_day == nday[imon]+1:
       if imon == len(months)-1:
          imon = 0
       else:
          imon += 1
       mon_day = 1
    mon = months[imon]
    day = mon_day

    #CREATE figures
    im1  = axarr[0,0].contourf(x,-1*zt,np.squeeze(btemp[i+toff,:,:]),  vt,   iedgecolors='face', extend = 'both')
    im2  = axarr[0,1].contourf(x,-1*zt,np.squeeze(mtemp[i+toff,:,:]),  vt,   iedgecolors='face', extend = 'both')
    im3  = axarr[1,0].contourf(x,-1*zvel,np.squeeze(bvel[i+toff,:,:]), vvel, cmap='bwr', iedgecolors='face', extend = 'both')
    im4  = axarr[1,1].contourf(x,-1*zvel,np.squeeze(mvel[i+toff,:,:]), vvel, cmap='bwr', iedgecolors='face', extend = 'both')
#    tx_str = str(year) + ' - ' + mon + ' - ' + str(day).zfill(2)
#    tx.set_text(tx_str)
    if i == 0:
#       axarr[0,0].invert_yaxis()
#       axarr[1,1].invert_yaxis() 
       plt.colorbar(im1, cax = cbar_ax_1,ticks = [tvals[0],np.mean(tvals),tvals[1]])
       plt.colorbar(im3, cax =cbar_ax_2, ticks = [vvals[0],np.mean(vvals),vvals[1]])
       cbar_ax_2.set_title('Horizontal Velocity (m/s)', rotation=270, va='center', ha='center', fontsize=16, y=.5)
       cbar_ax_1.set_title('Temperature ($^\circ$C)', rotation=270, va='center', ha='center', fontsize=16, y=.5)
       axarr[1,0].set_xlabel(xlabel[NB], x=1.1, fontsize=20,labelpad=30)
       axarr[1,0].set_ylabel('DEPTH (m)',y=1.1, fontsize=20,labelpad=30)
       axarr[0,0].set_title('HYCOM', fontsize=25, y = 1.02)
       axarr[0,1].set_title('MAPHIL MODEL', fontsize=25,y = 1.02)

    print i
#    print tx_str
##########
yr = 2009
imon = 10
mon_day = 3 

months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
nday   = [31,29,31,30,31,30,31,31,30,31,30,31]

grd = pyroms.grid.get_ROMS_grid('MaPhil')

rand_ncfile = '/t1/scratch/liz/tmpdir_MaPhil-LD.HCo06T/outputs/1996/MaPhil-LD.HCo06T_his_1996-02-01T12:00:00.nc'
bdfile = '/t1/scratch/liz/Inputs/MaPhil-LD.HCo05T/Boundary/hycom_BDRY.nc'
mafile = ['/t3/workdir/liz/MODELS/MAPHIL/PostProc/west_2009.nc',\
          '/t3/workdir/liz/MODELS/MAPHIL/PostProc/north_2009.nc',\
          '/t3/workdir/liz/MODELS/MAPHIL/PostProc/south_2009.nc']
#mafile = ['/t3/workdir/liz/MODELS/MAPHIL/PostProc/OBCFAC_edit/his/west_maphil.nc',\
#          '/t3/workdir/liz/MODELS/MAPHIL/PostProc/OBCFAC_edit/his/north_maphil.nc',\
#          '/t3/workdir/liz/MODELS/MAPHIL/PostProc/OBCFAC_edit/his/south_maphil.nc']




# FILE IDs
rfid = nc.Dataset(rand_ncfile)
bfid = nc.Dataset(bdfile)

# Variable names
btvar  = ['temp_west','temp_north','temp_south']
bvvar  = ['u_west','v_north','v_south']
mvvar  = ['u','v','v']

# TIME INDICES
#t1 = 1553 # 1/1/2014
t1 = 32 
tN = t1+60  # 
toff = 0 

# Template variables for creating z
rtemp = np.squeeze(rfid.variables['temp'][:])
time  = np.squeeze(bfid.variables['ocean_time'][t1:tN+1])
ndays = len(time)

# SLICE INDICES
ji = [0, grd.hgrid.y_v.shape[0]-1, 0]

# Contour limits
#tvals = [15,33]
tvals = [15,30]
vvals = [-1,1]
vt   = np.linspace(tvals[0], tvals[1], num=50,endpoint=True)
vvel = np.linspace(vvals[0], vvals[1], num=50,endpoint=True)

# Plot Labels
title_vals= ['WESTERN', 'NORTHERN', 'SOUTHERN']
xlabel = ['J POINTS', 'I POINTS', 'I POINTS']

# NB = [WEST, NORTH, SOUTH]
for NB in range(3):
 mfid = nc.Dataset(mafile[NB])

 # Template variables for creating z
 rvel  = np.squeeze(rfid.variables[mvvar[NB]][:])

 # Boundary Variables 
 btemp = np.squeeze(bfid.variables[btvar[NB]][t1:tN+1,:,:])
 bvel  = np.squeeze(bfid.variables[bvvar[NB]][t1:tN+1,:,:])

 # MAPHIL Values
 mtemp = np.squeeze(mfid.variables['temp'][:])
 mvel  = np.squeeze(mfid.variables[mvvar[NB]][:])
 if NB == 1:
    mvel = -1*mvel
    bvel = -1*bvel
 # Get Z values for depth
 if NB == 0: 
   slicet, zt, lon, lat = pyroms.tools.islice(rtemp, ji[NB], grd, Cpos='rho')
   slicev, zvel, lon, lat = pyroms.tools.islice(rvel,  ji[NB], grd, Cpos=mvvar[NB])
 else: 
   slicet, zt, lon, lat = pyroms.tools.jslice(rtemp, ji[NB], grd, Cpos='rho')
   slicev, zvel, lon, lat = pyroms.tools.jslice(rvel,  ji[NB], grd, Cpos=mvvar[NB])

 xmin,xmax = np.ma.notmasked_edges(slicet[0,:])

 # Plotting Along X
 #x   = np.tile(np.linspace(0, 1, num=zt[:].shape[1], endpoint=True),(zt[:].shape[0],1))
 x = np.tile(np.arange(zt.shape[1]),(zt[:].shape[0],1))

 # Create Figure
 f, axarr = plt.subplots(2,2, figsize=(14,14),sharex='col',sharey='row')
 for nt in range(2):
     axarr[nt,0].invert_yaxis()
     for ni in range(2):
         axarr[nt,ni].set_axis_bgcolor((0.7,0.7,0.7))
         axarr[nt,ni].set_xlim(xmin-1,xmax+1)
         axarr[nt,ni].set_ylim(340,0) 
 #CREATE cbar axes
 cbar_ax_1 = f.add_axes([0.93,0.55,0.03,0.33])
 cbar_ax_2 = f.add_axes([0.93,0.12,0.03,0.33])
 cbar_ax_2.text(.5, 1.1, 'IN', ha='center', va='center')
 cbar_ax_2.text(.5, -.1, 'OUT', ha='center', va='center')
 plt.subplots_adjust(hspace=0.1, wspace=0.1)
 title_t = title_vals[NB] + ' BOUNDARY' 
 f.suptitle(title_t, fontsize=30)
 # Run animation
 #ani = animation.FuncAnimation(f, updatefig,frames=np.len(time), blit=False)
 ani = animation.FuncAnimation(f, updatefig,frames=60, blit=False)
 str_nm = 'CT_rms_start_1_1_2014' + str(NB+1) + '.gif'
 ani.save(str_nm, writer = 'imagemagick',fps=4)
