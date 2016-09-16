import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

ncfile = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/1996/VIP-LD.HCo13T_avg_1996-05-17T00:00:00.nc'
fid = nc.Dataset(ncfile)
u = np.squeeze(fid.variables['u'][:])
v = np.squeeze(fid.variables['v'][:])
u_masked = np.ma.masked_where(u< 0.1, u)
u_masked.mask[:2,:,:]=True
u_masked.mask[:,:,:140]=True
u_masked.mask[:,:,440:]=True
u_masked.mask[:,:41,:]=True
u_masked.mask[:,290:,:]=True

for i in range(140,440): 
    for j in range(41,290):
        for k in range(1,u.shape[0]-1):
            if u_masked.mask[k,j,i] == False:
               t_val=0
               if u_masked[k,j,i+1]<0.1:    
                  t_val+=1
               if u_masked[k,j,i-1]<0.1:  
                  t_val+=1
               if u_masked[k,j+1,i]<0.1:      
                  t_val+=1
               if u_masked[k,j-1,i]<0.1:    
                  t_val+=1
               if u_masked[k+1,j,i]<0.1:      
                  t_val+=1
               if u_masked[k-1,j,i]<0.1:    
                  t_val+=1

               if (t_val<5):
                  u_masked.mask[k,j,i]=True
               if v[k,j,i]>0: 
                  u_masked.mask[k,j,i]=True 
plt.figure()
plt.pcolor(np.max(u_masked,axis=0))
plt.colorbar()
plt.show()
plt.xlim(0,u.shape[2])
plt.ylim(0,u.shape[1])
