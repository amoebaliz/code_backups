import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

#def get_curl(j,i):
#    ncfile = dir+tau_fil
#    fid = nc.Dataset(ncfile)
#    ws_curl = np.squeeze(fid.variables['ws_curl'][:,j,i])

#    return ws_curl

#def get_sst(j,i):
#    ncfile = dir+sst_fil
#    fid = nc.Dataset(ncfile)
#    sst = np.squeeze(fid.variables['sst'][:,j,i])
   
#    return sst

# GRID FILE SPECIFICS
VIP = pyroms.grid.get_ROMS_grid('VIP')

fid_1 = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/VIP-LD.HCo11T/TAU_Curl_96-99.nc')
fid_2 = nc.Dataset('/t3/workdir/liz/MODELS/VIP/PostProc/VIP-LD.HCo11T/sst_psi_96-99.nc')

R_val = np.zeros((341,601))

for j in range(341):
    print j
    for i in range(601):
        curl_ts = np.squeeze(fid_1.variables['ws_curl'][:,j,i])
        sst_ts  = np.squeeze(fid_2.variables['sst'][:,j,i])
        R = np.ma.corrcoef(curl_ts,sst_ts)
        R_val[j,i] = R[0,-1]

plt.figure()
plt.pcolor(R_val) 
plt.savefig('foo.png')
#plt.show()    

