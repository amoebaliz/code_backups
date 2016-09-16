import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as pltd

def convert_time(time):
    ref = dt.datetime(1900,1,1,0,0)
    date_vals = np.zeros(len(time))
    for nt in range(len(time)):
        day_time = ref + dt.timedelta(seconds=np.float(time[nt]))
        date_vals[nt] = pltd.date2num(day_time)
    return date_vals

def set_yaxes():

    # PADDING CALCULATIONS
    pltmin = np.mean(data_stor[0,:]) - stdfrac*np.std(data_stor[0,:])   
    submax = np.mean(data_stor[0,:]) + stdfrac*(2*(nplts-1)+1)*np.std(data_stor[0,:])

    pltmax = np.mean(data_stor[nplts-1,:]) + stdfrac*np.std(data_stor[nplts-1,:])
    submin = np.mean(data_stor[nplts-1,:]) - stdfrac*(2*(nplts-1)+1)*np.std(data_stor[nplts-1,:])

    lpd = (1.0*pltmin-np.min(data_stor[0,:]))/(submax-pltmin) + pdad
    upd = (1.0*np.max(data_stor[nplts-1,:])-pltmax)/(pltmax-submin) + pdad 

    # ITERATE OVER PLOTS    
    for nt in range(nplts):
        plt_rng = nplts*2*stdfrac*np.std(data_stor[nt,:]) 
        ymin = np.mean(data_stor[nt,:]) - stdfrac*(2*nt+1)*np.std(data_stor[nt,:])
        ymax = np.mean(data_stor[nt,:]) + stdfrac*(2*(nplts - (nt+1))+1)*np.std(data_stor[nt,:])

        # APPLY PADS
        ymin = ymin - lpd*plt_rng
        ymax = ymax + upd*plt_rng

        # LABEL POSITION CENTERED AT 0.0 TICK 
        ypos = (0.0 - ymin)/(ymax-ymin)

        axs[nt].set_ylim(ymin,ymax)
        axs[nt].set_yticks([-1*tcks[nt],0,tcks[nt]])
        axs[nt].set_ylabel(ylabels[nt], y = ypos, rotation = 90 + (nt%2)*180, fontsize=16, labelpad = 15)     
        if nt%2==1:
           axs[nt].yaxis.tick_right()
           axs[nt].yaxis.set_label_position("right")
       
        else:
           axs[nt].yaxis.tick_left()
           axs[nt].yaxis.set_label_position("left")

    #ax.figure.canvas.draw() 
# ------------------------------------------------------ #
# DATA DETAILS
dir = '/t3/workdir/liz/MODELS/VIP/PostProc/'
fils = ['vip_transports/vip_u_100m_transport_col278_rows144-178.nc',\
        'zeta_278_96_99.nc',\
        'SSH_1996-1999.nc', \
        'sustr_278_96_99.nc']

vars = ['_transport','zeta','zeta','sustr']

# CONSTANTS
fSv = 1000000

# PLOT SPECS
nplts = len(fils)
colval = ['-b','-g','-k','-r']
ylabels = []
zcent = [1,1,1]
tcks = [.8,.04,4,0.01]
xmin  = pltd.date2num(dt.datetime(1996,1,1,0,0))
xmax  = pltd.date2num(dt.datetime(1999,12,31,0,0))
ylabels = ['Transport (Sv)',r'$\Delta$ SSH (cm)',r'$\Delta$ SSH (cm)','thing2 ()']
stdfrac = 3 # # OF STDEVs (*.5) DEVOTED TO PLOT SPACE
pdad = .01


fig, ax = plt.subplots(1, figsize=(10,10))
[fig.add_axes(ax.twinx()) for n in range(nplts-1)]
axs = fig.get_axes()
# PLOTTING and STORING
for nt in range(nplts+1):
    ncfile = dir + fils[nt]
    fid = nc.Dataset(ncfile)

    if nt == 0: 
       # Get time for shared x-axis
       time = fid.variables['ocean_time'][:]
       date_vals = convert_time(time)
       
       # TRANSPORT variable modifications
       var =  (np.squeeze(fid.variables['pos_transport'][:]) + \
               np.squeeze(fid.variables['neg_transport'][:]))/fSv
       data_stor = np.empty([nplts, len(date_vals)])
    else: 
       var = np.squeeze(fid.variables[vars[nt]][:])
       if nt == 1:
          # SSH variable modification
          #var = 100*(np.squeeze(np.mean(var[:,65:149,180],axis=1))-np.squeeze(np.mean(var[:,149:199,318],axis=1)))
          var = 100*(np.squeeze(np.mean(var[:,:,0],axis=1))-np.squeeze(np.mean(var[:,:,1],axis=1)))
       elif nt == 2: 
          print var.shape
          var = 100*(np.squeeze(np.mean(var[:,65:149,180],axis=1))-np.squeeze(np.mean(var[:,149:199,318],axis=1)))
       else:
          # TauX variable modification
          var = np.mean(np.mean(var,axis=2),axis=1)

    axs[nt].plot(date_vals,var,colval[nt])
    data_stor[nt,:] = var 

data_stor = data_stor.reshape(nplts,-1)
set_yaxes()

#print np.max(np.correlate(data_stor[0,:],data_stor[1,:],mode='full'))
print 'local:', np.corrcoef(data_stor[0,:],data_stor[1,:])[0,1]
print 'distal:',np.corrcoef(data_stor[0,:],data_stor[2,:])[0,1]
#       print np.corrcoef(transp[:365],dvar[:365])[0,1]
#       print np.corrcoef(transp[1:365],dvar[:365-1])[0,1]
#       print np.corrcoef(transp[2:365],dvar[:365-2])[0,1]
#       print np.corrcoef(transp[:365-1],dvar[1:365])[0,1]
#       print np.corrcoef(transp[:365-2],dvar[2:365])[0,1]
#auto_cor = np.correlate(transp-np.mean(transp),dvar-np.mean(dvar),mode='full')
#      print len(transp)-np.where(auto_cor == np.min(auto_cor))[0][0]


#       print transp.shape
#       print dvar.shape
# print np.corrcoef(transp,transp1)[0,1]

# TEXT LABELS
#txpos = pltd.date2num(dt.datetime(1900,1,1,0,0))-30
#ax.text(txpos, 30.2, 'CoRTAD Climatology', color = 'g', fontsize=16)
#ax.errorbar(txpos, 29.2, yerr = .5, fmt='ok', ms = 2, ecolor=[.5,.5,.5],elinewidth=2,capsize=0)
#ax.text(txpos+5, 29.2, '1998', va = 'center', color = 'k', fontsize=16)

#ax2.text(txpos, Y1[0]+2, 'MERRA Climatology', fontsize=16)
#ax2.text(txpos, Y2[0], '1998 Anomaly', color = 'b', fontsize=16)

# X AXIS
ax.xaxis_date()
ax.set_xlim(xmin,xmax)
ax.xaxis.set_ticks_position('bottom')
# TEXT LABELS

# FINAL DETAILS
plt.tight_layout(pad=2)
plt.savefig('drivers.png')
plt.show()
