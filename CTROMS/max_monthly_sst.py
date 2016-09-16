import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

dir = '/t3/workdir/liz/MODELS/CTROMS/PostProc/'

frame_len = 31
years = range(2000,2005)
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]

ilats = range(400,450)
ilons = range(400,450)

stored_run_max = np.zeros((len(years),len(ilats),len(ilons)))
stored_mon_max = np.zeros((len(years),len(ilats),len(ilons)))
for yr in range(len(years)):

    print years[yr]
    # Load 1yr SST data
    if years[yr]%4==0:
       ndays[1] = 29
    else:
       ndays[1] = 28

    file = dir + str(years[yr]) + '_CT_SST.nc'
    fid = nc.Dataset(file)
    temp = np.squeeze(fid.variables['temp'][:,:,ilats,ilons])
    # CREATE and PRINT MONTHLY MEANS
    mon_max = np.zeros((12,len(ilats),len(ilons)))
    iday=0 
    for nmon in range(12):
        mon_max[nmon,:,:] = np.mean(temp[iday:(iday+ndays[nmon]),:,:],axis=0) 
        iday+=ndays[nmon]

    stored_mon_max[yr,:,:] = np.max(mon_max,axis=0)

    # EVALUATE 31 DAY PERIODS WITHIN 1YR; STORE MAX
    for nframe in range(temp.shape[0]-frame_len+1):
        run_max = np.mean(temp[nframe:nframe+frame_len,:,:],axis = 0)
        stored_run_max[yr,:,:] = np.where(run_max>stored_run_max[yr,:,:],run_max,stored_run_max[yr,:,:])


print 'MAX MEAN MON'
print np.mean(stored_mon_max[:,0,0],axis=0)

print 'RUNNING 31-DAY MEAN MAX'
print np.mean(stored_run_max[:,0,0],axis=0)



