import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as pltd
import datetime as dt

# SSH
ssh_file ='/t3/workdir/liz/MODELS/VIP/PostProc/VIP13_SSH_1996-1999.nc'
ssh_fid = nc.Dataset(ssh_file)
ssh_A = np.mean(np.squeeze(ssh_fid.variables['zeta'][:,65:149,180]),axis=1)
ssh_B = np.nanmean(np.squeeze(ssh_fid.variables['zeta'][:,146:199,318]),axis=1)


#A = [172,116];    store_A = np.squeeze(ssh[:,A[1],A[0]])
#B = [175,75];     store_B = np.squeeze(ssh[:,B[1],B[0]])
#C = [224,127];    store_C = np.squeeze(ssh[:,C[1],C[0]])
#D = [265,154];    store_D = np.squeeze(ssh[:,D[1],D[0]])
#E = [323,175];    store_E = np.squeeze(ssh[:,E[1],E[0]])
#F = [383,202];    store_F = np.squeeze(ssh[:,F[1],F[0]])

#ssh_avg = np.mean(ssh,axis=0)
ssh_dif = ssh_A-ssh_B
print ssh_B
# TRANSPORT
# TRANSECT POINTS [i,j]
vip1 = [278,144]; vip2 = [278,178]
sv_file ='vip_u_ALLm_transport_col278_rows144-178.nc'
sv_fid = nc.Dataset(sv_file)
pos = sv_fid.variables['pos_transport'][:]
neg = sv_fid.variables['neg_transport'][:]
fSv = 1000000
total_trans = (pos+neg)/fSv

T =  np.corrcoef(total_trans, ssh_dif)
print T

# TIME
time = ssh_fid.variables['ocean_time'][:]
ref = dt.datetime(1900,1,1,0,0)
date2 = np.zeros(len(time))
for nt in range(len(time)):
    date_val = ref + dt.timedelta(seconds=time[nt])
    date2[nt] = pltd.date2num(date_val)

# PLOT MAP
#f = plt.figure()
#ax = f.add_subplot(111,axisbg=[0.7,0.7,0.7])
#P = ax.pcolor(ssh_avg,vmin=0.72,vmax=.75)
#plt.colorbar(P)
#plt.plot(A[0],A[1],'o'); plt.text(A[0]+3,A[1]+3,'A')
#ax.plot(B[0],B[1],'o'); plt.text(B[0]+3,B[1]+3,'B')
#plt.plot(C[0],C[1],'o'); plt.text(C[0]+3,C[1]+3,'C')
#plt.plot(D[0],D[1],'o'); plt.text(D[0]+3,D[1]+3,'D')
#plt.plot(E[0],E[1],'o'); plt.text(E[0]+3,E[1]+3,'E')
#ax.plot(F[0],F[1],'o'); plt.text(F[0]+3,F[1]+3,'F')
#ax.plot([vip1[0],vip2[0]],[vip1[1],vip2[1]],'-k',linewidth=2.0)
#ax.set_xlim([0,ssh_avg.shape[1]])
#ax.set_ylim([0,ssh_avg.shape[0]])

#plt.figure()
#A, = plt.plot(nums,store_A,'-k')
#B, = plt.plot(nums,store_B,'-b')
#C, = plt.plot(nums,store_C,'-g')
#D, = plt.plot(nums,store_D,'-r')
#E, = plt.plot(nums,store_E,'-y')
#F, = plt.plot(nums,store_F,'-r')
#plt.legend([A,B,C,D,E,F],['A','B','C','D','E','F'])

# PLOT SSH DIF AND TRANSPORT

fig = plt.figure(figsize=(15, 3)) 
ax1 = fig.add_subplot(111)
dif_plot, = ax1.plot(date2,ssh_dif,'-k',linewidth=1.5)
ax1.set_ylim([-.06,.06])
ax1.set_yticks([-.06,0,.06])
ax1.set_ylabel('SSH difference (m)')
#ax1.set_xlabel('Time (days)')
ax1.set_xlim([date2[0],date2[-1]])
ax1.xaxis_date()

ax2 = ax1.twinx()
tot_trans, = ax2.plot(date2,total_trans,'-b')
ax2.plot([date2[0],date2[-1]],[0,0],'--k')
ax2.set_ylim([-.7,.7])
ax2.set_yticks([-.7,0,.7])
ax2.set_ylabel('Transport (Sv)',color = 'b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')

#ax2.text(date2[90],.5,'R = %.2f' % T[-1,0])
#plt.legend([dif_plot,tot_trans],['B-F','transport'],loc = 3)
#plt.title('Comparison of SSH difference between points B and F \n and zonal surface transport through the VIP')
plt.savefig('b.png')
#plt.show()
