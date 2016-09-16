import netCDF4 as nc
import datetime as dt
import numpy as np
########################## INITIAL CONDITION ######################
# NOTE: copied original anv file created by make_ini_file and renamed it
fid = nc.Dataset('./VIP_Ini_fromCORAL_y2007m01md02.nc','a')
#fid = nc.Dataset('./MaPhil_Ini_fromCORAL_y1978m02md18.nc','a')
time = fid.variables['ocean_time']
tmp = time[:]

ref = dt.datetime(1900,1,1,0,0)
newdate = ref + dt.timedelta(seconds=tmp[0])
delta = newdate - ref
toto = delta.days + (delta.seconds / 86400.)

time.units = 'days since 1900-01-01 00:00:00'
time[:] = toto
fid.close()
############### BOUNDARY CONDITION ####################
fid = nc.Dataset('./VIP_BRY_y2007.nc','a') # read mode
#fid = nc.Dataset('./MaPhil_BRY_y1978.nc','a')
time = fid.variables['ocean_time']
tmp = time[:]
timeout = tmp.copy()

ref = dt.datetime(1900,1,1,0,0)

for kt in np.arange(len(tmp)):
newdate = ref + dt.timedelta(seconds=tmp[kt])
delta = newdate - ref
toto = delta.days + (delta.seconds / 86400.)
timeout[kt] = delta.days + (delta.seconds / 86400.)


print timeout

time.units = 'days since 1900-01-01 00:00:00'
time[:] = timeout
fid.close()

############### NUDGING CONDITION ####################
fid = nc.Dataset('./coral_clim_2007_VIP.nc','a')
#fid = nc.Dataset('./MaPhil_BRY_y1978.nc','a')
time = fid.variables['ocean_time']
tmp = time[:]
timeout = tmp.copy()

ref = dt.datetime(1900,1,1,0,0)

for kt in np.arange(len(tmp)):
newdate = ref + dt.timedelta(seconds=tmp[kt])
delta = newdate - ref
toto = delta.days + (delta.seconds / 86400.)
timeout[kt] = delta.days + (delta.seconds / 86400.)


print timeout

time.units = 'days since 1900-01-01 00:00:00'
time[:] = timeout
fid.close()

########



fid = nc.Dataset('./coral_monthly_avg_2008_01_VIP.nc','a')
#fid = nc.Dataset('./MaPhil_BRY_y1978.nc','a')
time = fid.variables['ocean_time']

ref = dt.datetime(1900,1,1,0,0)
newdate = dt.datetime(2008,1,16,12,0)
delta = newdate - ref
timeout = delta.days*(24*60*60)+delta.seconds

time.units = 'seconds since 1900-01-01 00:00:00'
time[:] = timeout
fid.close()



# IC jan 02 1979
	#./coral_avg_07673.nc
# Boundary files start at: ./coral_avg_07673.nc
		 # end at: ./coral_avg_07703.nc 

######## CONCATINATE BOUNDARY CONDITIONS ##########
ncrcat *bdry.nc -o VIP_BRY_y1979.nc # fro; jan 01 1979
ncrcat *bdry.nc -o VIP_BRY_y2007.nc 
ncrcat *bdry.nc -o MaPhil_BRY_y1978.nc


############ MAKE INITIAL CONDITION COPY
cp coral_avg_07674_VIP.nc VIP_Ini_fromCORAL_y1979m01md02.nc
cp coral_avg_17901_VIP.nc VIP_Ini_fromCORAL_y2007m01md02.nc
cp coral_avg_07356_MaPhil.nc MaPhil_Ini_fromCORAL_y1978m02md18.nc


####### Determine Ocean_time

#fid = nc.Dataset('/data/external/P2/ROMS/CORAL/RUN16/coral_avg_17901.nc','r')
fid = nc.Dataset('./coral_monthly_avg_2008_01_VIP.nc','r')
time = fid.variables['ocean_time']
tmp = time[:]
ref = dt.datetime(1900,1,1,0,0)
newdate = ref + dt.timedelta(seconds=tmp[0])
print newdate



for year in np.arange(2010,2010+1):
YEAR = str(year).zfill(4)
for month in np.arange(1,12+1):
MONTH = str(month).zfill(2)
ncfile = 'hycom_' + YEAR + '_' + MONTH + '_clim_MaPhil.nc'

ncfile ='hycom_clim.nc'
fid = nc.Dataset(ncfile,'r')
time = fid.variables['ocean_time']
tmp = time[:]
#print tmp[0]
ref = dt.datetime(1900,1,1,0,0)
for i in np.arange(len(tmp)):
newdate = ref + dt.timedelta(days=tmp[i])
print newdate


#######################
import netCDF4 as nc
import numpy as np
import datetime as dt

ref = dt.datetime(1900,1,1,0,0)
nums = range(17900,18262+1)
timeout = np.zeros(len(nums))
for nt in nums:
ncfile = '/data/external/P2/ROMS/CORAL/RUN16/coral_avg_'+str(nt)+'.nc'
fid = nc.Dataset(ncfile,'r')
ocean_time= fid.variables['ocean_time'][:]
newdate = ref + dt.timedelta(seconds=ocean_time[0])
delta = newdate - ref
toto = delta.days + (delta.seconds / 86400.)
timeout[nt-nums[0]]=toto

print timeout



