# calculate CTROMS transports through VIP
import netCDF4 as nc
import numpy as np
import datetime as dt
import pyroms as py
import pyroms.tools as pyt

cols = np.array([192,240,280,341])
rows = np.array([[70,148],[95,152],[144,178],[140,200]])
h1=0; h2=-50

# TIME RANGE 1996-1999
#ct_nums = range(13881,15342+1)

vip_dir = '/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo11T/'
VIP = py.grid.get_ROMS_grid('VIP')

def make_ncfile(i):
    file = 'vip_u_transport_col' + str(cols[i]) + '_96-99.nc'
    ncid = nc.Dataset(file,'w')
    ncid.createDimension('ocean_time', None)
    ncid.title = 'zonal transports between ' + str(h1) + ':' + str(-1*h2) + ' meters, at vip column '+ str(cols[i]) + ' and rows ' + str(rows[i,0]) + ':' + str(rows[i,1])

    ocean_time = ncid.createVariable('ocean_time', 'f8',('ocean_time',))
    ocean_time.units = 'seconds since 1900-01-01 00:00:00'

    pos_transport = ncid.createVariable('pos_transport','f8',('ocean_time',))
    pos_transport.units = 'm3/s'
    pos_transport.long_name = 'eastward transport'

    neg_transport = ncid.createVariable('neg_transport','f8',('ocean_time',))
    neg_transport.units = 'm3/s'
    neg_transport.long_name = 'westward transport'

    ncid.close()
    return file

def get_transports(i,j,yr,mon,day,file):
    yrstr = str(yr)
    monstr = str(mon).zfill(2)
    daystr = str(day).zfill(2) 

    ncfile = vip_dir + yrstr + '/VIP-LD.HCo11T_avg_' + yrstr + '-' + monstr + '-' + daystr + 'T00:00:00.nc'
    fid = nc.Dataset(ncfile,'r')

    u = fid.variables['u'][:]; u.setflags(write=False)
    v = fid.variables['v'][:]
    time = fid.variables['ocean_time'][:]

    v = np.squeeze(v)
    posu = np.squeeze(np.ma.copy(u))
    negu = np.squeeze(np.ma.copy(u))
    posu[posu<0]=0
    negu[negu>0]=0

    pos_ut,transpv = pyt.section_transport_z(posu, v, VIP, cols[j], cols[j], rows[j,0], rows[j,1], h1, h2)
    neg_ut,transpv = pyt.section_transport_z(negu, v, VIP, cols[j], cols[j], rows[j,0], rows[j,1], h1, h2)

    ncid = nc.Dataset(file, 'a')
    ocean_time = ncid.variables['ocean_time'];       ocean_time[i] = time
    pos_transport = ncid.variables['pos_transport']; pos_transport[i] = pos_ut
    neg_transport = ncid.variables['neg_transport']; neg_transport[i] = neg_ut
    ncid.close()

ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
for j in range(len(cols)):
    file = make_ncfile(j)
    n = 0     
    for yr in range(1996,1999+1):

        if (yr % 4) == 0:
            ndays[1]=29
            print yr, 'LEAP'
        else: 
            ndays[1]=28
            print yr

        for mon in range(1,12+1):
            for day in range(1,ndays[mon-1]):
                get_transports(n,j,yr,mon,day,file) 
                 
                n+=1
     
