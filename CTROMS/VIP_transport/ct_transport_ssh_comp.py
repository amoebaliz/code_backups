import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt
import pyroms as py
import pyroms.tools as pyt

SSH_DIR = '/t3/workdir/liz/ANALYSES/CTROMS/SSH_VIP/'
TRANS_DIR = '/t3/workdir/liz/ANALYSES/CTROMS/TRANSPORT_VIP/'

#nc1 = SSH_DIR + 'ct_zeta_row476_col360.nc'
#nc1 = SSH_DIR + 'ct_zeta_row468_col378.nc'
#nc1 = SSH_DIR + 'ct_zeta_row469_col378.nc'
#nc1 = SSH_DIR + 'ct_zeta_row474_col379.nc'
#nc1 = SSH_DIR + 'ct_zeta_row471_col382.nc'
nc1 = SSH_DIR + 'ct_zeta_row471_col385.nc'
#nc1 = SSH_DIR + 'ct_zeta_row470_col390.nc'


#nc2 = SSH_DIR + 'ct_zeta_row471_col393.nc'
nc2 = SSH_DIR + 'ct_zeta_row471_col397.nc'
#nc2 = SSH_DIR + 'ct_zeta_row471_col401.nc'
#nc2 = SSH_DIR + 'ct_zeta_row472_col402.nc'
#nc2 = SSH_DIR + 'ct_zeta_row464_col430.nc'

#nc3 = TRANS_DIR + 'ctroms_u_transport_col395.nc'
nc3 = TRANS_DIR + 'ctroms_u_transport_col392.nc'
#nc3 = TRANS_DIR + 'ctroms_u_transport_col384.nc'

#nc4 = 'ct_roms_avg_2007_zeta.nc'
#nc4 = '/data/external/P2/ROMS/CORAL/RUN16/coral_avg_11054.nc'

fid1 = nc.Dataset(nc1,'r')
fid2 = nc.Dataset(nc2,'r')
fid3 = nc.Dataset(nc3,'r')
#fid4 = nc.Dataset(nc4,'r')


zeta1 = np.squeeze(fid1.variables['zeta'][732:])
zeta2 = np.squeeze(fid2.variables['zeta'][732:])
#zeta_mean = np.squeeze(fid4.variables['zeta'][:])

pos_trans = fid3.variables['pos_transport'][732:]
neg_trans = fid3.variables['neg_transport'][732:]
date_val = fid3.variables['ocean_time'][732:]

ssh_diff = zeta1-zeta2

fSv = 1000000
sum_trans = (pos_trans+neg_trans)/fSv

ref = dt.datetime(1900,1,1,0,0)
date = np.zeros(len(ssh_diff))
for nt in range(len(ssh_diff)):
   day_time = ref + dt.timedelta(seconds=date_val[nt])
   date[nt] = pltd.date2num(day_time)


T = np.corrcoef(ssh_diff,sum_trans)
print T
fig,ax1 = plt.subplots(1,figsize = (16,5))

dif_plot, = ax1.plot(date,ssh_diff,color = [0.2,0.2,0.2],linewidth=2.0)
ax1.plot([date[0],date[-1]],[0,0],'--k')
ax1.set_xlim(date[0],date[-1])
ax1.set_ylim([-.15,.15])
ax1.set_ylabel('SSH difference (m)')
ax1.set_xlabel('Time (days)')
ax1.xaxis_date()
ax1.set_title('Comparison of CT-ROMS SSH difference and zonal surface transport in the VIP')

ax2 = ax1.twinx()
tot_trans, = ax2.plot(date,sum_trans,'-r')
ax2.set_ylabel('Transport (Sv)',color = 'r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

ax1.text(date[100],.1,'R = %.2f' % T[-1,0])
plt.legend([dif_plot,tot_trans],['A-B','transport'],loc = 1)

#plt.figure()

#plt.pcolor(zeta_mean)
#plt.colorbar()
#plt.plot(float(nc1[18:21]),float(nc1[11:14]),'o') 
#plt.text(float(nc1[18:21])+2,float(nc1[11:14])+2,'A')

#plt.plot(float(nc2[18:21]),float(nc2[11:14]),'o') 
#plt.text(float(nc2[18:21])+2,float(nc2[11:14])+2,'B')

#plt.plot([395,395],[469,474],'-k',linewidth=2.0)
#plt.clim([0.74, 0.78])
#plt.xlim([300,450])
#plt.ylim([400,640])
plt.show()


