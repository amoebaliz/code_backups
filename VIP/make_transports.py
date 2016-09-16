import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as pltd
import datetime as dt
import pyroms as py
import pyroms.tools as pyt

ref = dt.datetime(1900,1,1,0,0)

# TRANSECT POINTS [i,j]
vip1 = [278,144]; vip2 = [278,178]
#vip1 = [280,145]; vip2 = [280,178]
#vip1 = [318,146]; vip2 = [318,199]
# DEPTH RANGES
h1 = 0;  h2 = -50

VIP = py.grid.get_ROMS_grid('VIP')
vip_dir = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/'
vip_nums = 365*4+1

def get_transports(ncfile):
    fid = nc.Dataset(ncfile,'r')
    u = fid.variables['u'][:]; u.setflags(write=False)
    v = fid.variables['v'][:]        
    time = fid.variables['ocean_time'][:]
    
    v = np.squeeze(v)
    posu = np.squeeze(np.ma.copy(u)) 
    negu = np.squeeze(np.ma.copy(u))   
    posu[posu<0]=0 
    negu[negu>0]=0
    day_time = ref + dt.timedelta(seconds=time[0])
    
    pos_ut,transpv = pyt.section_transport_z(posu, v, VIP, vip1[0], vip2[0], vip1[1], vip2[1], h1, h2)
    neg_ut,transpv = pyt.section_transport_z(negu, v, VIP, vip1[0], vip2[0], vip1[1], vip2[1], h1, h2)
    date2 = pltd.date2num(day_time)
    return pos_ut, neg_ut, time, date2

# Building Plot values 
pos_u_trans = np.zeros(vip_nums)
neg_u_trans = np.zeros(vip_nums)
date2 = np.zeros(vip_nums)
delta_time = np.zeros(vip_nums)

ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
nt=0
for nyr in range(1996,1999+1):
    ineq = nyr%4
    if ineq == 0:
       ndays[1] = 29   
    else:
       ndays[1] = 28
    for nmon in range(1,12+1):
        for nday in range(1,ndays[nmon-1]+1):
            file = vip_dir + str(nyr) + '/VIP-LD.HCo13T_avg_' + str(nyr)+'-'+str(nmon).zfill(2)+'-'+str(nday).zfill(2)+'T00:00:00.nc'
            pos_u_trans[nt], neg_u_trans[nt], delta_time[nt], date2[nt]  = get_transports(file)[:]
            nt+=1
print 'save to file'

# SAVE TO NETCDF file
file = 'vip_u_50m_transport_col'+str(vip1[0])+'_rows'+str(vip1[1])+'-'+str(vip2[1])+'.nc'
ncid = nc.Dataset(file,'w')
ncid.createDimension('ocean_time', None)
ncid.title = 'zonal transports between ' + str(h1) + ':' + str(-1*h2) + ' meters, at VIP column '+str(vip1[0])+' and rows '+str(vip1[1])+':'+str(vip2[1])

ocean_time = ncid.createVariable('ocean_time', 'f8',('ocean_time',))
ocean_time.units = 'seconds since 1900-01-01 00:00:00'
ocean_time[:] = delta_time[:]

pos_transport = ncid.createVariable('pos_transport','f8',('ocean_time',))
pos_transport.units = 'm3/s'
pos_transport.long_name = 'eastward transport'
pos_transport[:] = pos_u_trans[:]

neg_transport = ncid.createVariable('neg_transport','f8',('ocean_time',))
neg_transport.units = 'm3/s'
neg_transport.long_name = 'westward transport'
neg_transport[:] = neg_u_trans[:]
ncid.close()
