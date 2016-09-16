import numpy as np
import netCDF4 as nc

ct_files = ['ct_VIPbound_sst_WEEKLY_ts_A.nc','ct_VIPbound_sst_WEEKLY_ts_B.nc','ct_VIPbound_sst_WEEKLY_ts_C1.nc','ct_VIPbound_sst_WEEKLY_ts_C2.nc','ct_VIPbound_sst_WEEKLY_ts_D.nc']

cortad_files = ['cortad_VIPbound_sst_ts_A.nc','cortad_VIPbound_sst_ts_B.nc','cortad_VIPbound_sst_ts_C1.nc','cortad_VIPbound_sst_ts_C2.nc','cortad_VIPbound_sst_ts_D.nc']

for nt in range(len(ct_files)):

    fid1 = nc.Dataset(ct_files[nt])
    fid2 = nc.Dataset(cortad_files[nt])
    ct_sst = np.squeeze(fid1.variables['temp'][:])
    cortad_sst = np.squeeze(fid2.variables['sst'][:])
    
    avg_dif = np.mean(ct_sst-cortad_sst)
    bound_std = np.std(ct_sst-cortad_sst)

    print avg_dif
    print bound_std
