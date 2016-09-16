import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

## OBJECTIVES:
#     3 part plots of winds, ssh, VIP velocity/transport
#     1st: plot all years faded; mean in darker color (more compact)
#     2nd: plot entire time series (can indicate enso years...)

def plot_winds():

def plot_ssh():

def plot_current():
    

# FILES
wind_file = '' 
wind_fid = nc.Dataset(wind_file)

ssh_file = ''
ssh_fid = nc.Dataset(ssh_file)

current_file = '' 
current_fid = nc.Dataset(current_file)

soi_file = ''

# EXTRACT VARIABLES
= wind_fid.variables[''][:]
= wind_fid.variables[''][:]
= ssh_fid.variables[''][:]
= current_fid.variables[''][:]

# ENTIRE TIME SERIES
plt.figure()
# plot winds, ssh, current, soi (?)

# CONDENSED TIME SERIES 
fig2 = plt.figure()
ax2 = 

windu_stored = np.zeros((len(yrs),365))
windv_stored = np.zeros((len(yrs),365))
ssh_stored = np.zeros((len(yrs),365))
current_stored = np.zeros((len(yrs),365))

for nyr in range(len(yrs)):
    windu_stored[nyr,:],windv_stored[nyr,:]= plot_winds(nyr)
    ssh_stored[nyr,:]= plot_ssh(nyr)
    current_stored[nyr,:]= plot_current(nyr)

# plot mean
ax2[0].quiver
ax2[1].plot(np.mean(ssh_stored,axis = 0))
ax2[2].plot(np.mean(current_stored,axis = 0))
 
plt.show()
