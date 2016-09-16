import numpy as np
import netCDF4 as nc
from tempfile import TemporaryFile


def get_manila_tides(i):
    man_tids = []
    for nt in range(i):
        tid_fil = tid_dir + tidal_files[nt]
#        if nt == 0:
#           arr_tmp=np.loadtxt(tid_fil,skiprows=489,usecols=(3,4,5,6,7,8,9,10,11,12,13,14))
#        else:
        arr_tmp=np.loadtxt(tid_fil,skiprows=1,usecols=(3,4,5,6,7,8,9,10,11,12,13,14))
        arr = arr_tmp.flatten()
        arr = np.ma.masked_where(arr > 9000,arr)
        man_tids = np.ma.append(man_tids,arr)
    tid_mean = np.mean(man_tids)
    arr_norm_m = (man_tids - tid_mean)/1000.

    return arr_norm_m

a = 11568
b = -820
#a =0
#b = -1
sta_fil ='/t1/scratch/liz/tmpdir_VIP-LD.HCo13T/VIP-LD.HCo13T_sta.nc'
fid = nc.Dataset(sta_fil)
vip_ssh = fid.variables['zeta'][a:b:4,14]
vip_ssh = vip_ssh - np.mean(vip_ssh)

tidal_files = ['h370a96.dat','h370a97.dat','h370a98.dat','h370a99.dat']
tid_dir =  '/data/external/P1/Data/tide_gauge/Manila/'
arr_norm_m = get_manila_tides(len(tidal_files))

print len(vip_ssh)/24

R_store = np.empty(len(vip_ssh)/24)
R_store.fill(np.nan)
outfile = TemporaryFile()
for nt in range(len(R_store)):
    c = 24*nt
    d = c + 30*24

    if ((c-30*24>= 0) and (d<=len(R_store)*24)):
    
       vssh = vip_ssh[c:d]
       mssh = arr_norm_m[c:d]
       R_val = np.ma.corrcoef(vssh,mssh)
       R_store[nt] = R_val[0,-1]

np.savetxt('daily_tidal_R_values',R_store)
