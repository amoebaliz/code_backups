import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dt

def get_taux(j,i):
    taux = []
    for n_yr in range(98,99+1):
#    for n_yr in range(96,99+1):
        file = dir + 'TAUX_' + str(n_yr) + '.nc'
        fid  = nc.Dataset(file)
        tau1 = np.squeeze(fid.variables['sustr'][:,j+1,i])
        tau2 = np.squeeze(fid.variables['sustr'][:,j+1,i+1])
        tau  = np.mean((tau1,tau2),axis = 0)
        taux  = np.concatenate([taux,tau])      
    return taux

def get_ssh(j,i):
    ssh = []
    for n_yr in range(98,99+1):
#    for n_yr in range(96,99+1):
        file = dir + 'SSH_' + str(n_yr) + '.nc'
        fid  = nc.Dataset(file)
        zeta = np.squeeze(fid.variables['zeta'][:,j+1,i+1])
        ssh  = np.concatenate([ssh,zeta]) 
     
    return ssh

###

dir = '/t3/workdir/liz/MODELS/VIP/PostProc/tmpdir_VIP-LD.HCo11T/' 

####### Generate correlation grid
plot_corr = np.zeros((340,600))

for i in range(plot_corr.shape[1]):
    for j in range(plot_corr.shape[0]):
        ssh = get_ssh(j,i)
        taux = get_taux(j,i)
        plot_corr[j,i] = np.corrcoef(taux,ssh)[0,1]

####### Generate Figure 

fig = plt.figure()
plt.pcolor(plot_corr)
plt.show()


