import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

######################

ct_trans_file = '/t3/workdir/liz/MODELS/CTROMS/TRANSPORT_VIP/ctroms_u_transport_col384.nc'
fid = nc.Dataset(ct_trans_file)
pos_trans = np.squeeze(fid.variables['pos_transport'][732:])
neg_trans = np.squeeze(fid.variables['neg_transport'][732:])
kSv = 1e6
yrs = np.arange(1960,2007+1) 

#####################

all_pos_subs = np.empty((len(yrs),365))
all_pos_subs[:] = np.NAN
all_neg_subs = np.empty((len(yrs),365))
all_neg_subs[:] = np.NAN

istart=0

plt.figure()
for nt in yrs:

    # Leap years
    if yrs[nt-yrs[0]]%4==0:
       ndays = 366
    else:
       ndays = 365

    # set end bounds
    if nt == 2007:
       iend = None

    else:
       iend=istart+ndays
    # define and store annual data
    pos_plot = pos_trans[istart:iend]
    neg_plot = neg_trans[istart:iend]
    if len(pos_plot)>365:
        pos_plot = pos_plot[:365]
        neg_plot = neg_plot[:365]

    all_pos_subs[nt-yrs,:len(pos_plot)] = pos_plot
    all_neg_subs[nt-yrs,:len(neg_plot)] = neg_plot
    
    # plot
#    if yrs[nt-yrs[0]] == 1999:
#       pos_cval = [.2,.2,.2]
#       neg_cval = [.2,.2,.2]
#    else: 
    pos_cval = [.9,.6,.6]
    neg_cval = [.6,.6,.9]
    plt.plot(np.array(range(len(pos_plot))),pos_plot/kSv, c = pos_cval)
    plt.plot(np.array(range(len(neg_plot))),neg_plot/kSv, c = neg_cval)  

    # update start index
    istart=iend
plt.plot(np.array(range(365)),np.mean(all_pos_subs,axis=0)/kSv, c = [1,0,0])
plt.plot(np.array(range(365)),np.mean(all_neg_subs,axis=0)/kSv, c = [0,0,1])
plt.xlim(0,366)
plt.show()










