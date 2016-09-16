import netCDF4 as nc
import numpy as npy
import matplotlib.pylab as plt
from matplotlib import cm as cm

class CCS_index():

        def __init__(self,ncgrid):
                self.ncgrid = ncgrid
                self.read_grid()


        def read_grid(self):
                fid = nc.Dataset(self.ncgrid,'r')
                self.lonrho = fid.variables['lon_rho'][:]
                self.latrho = fid.variables['lat_rho'][:]
                self.lonvert = fid.variables['lon_vert'][:]
                self.latvert = fid.variables['lat_vert'][:]
                self.maskrho = fid.variables['mask_rho'][:]
                fid.close()
                return None

        def distance_on_unit_sphere(self,lat1, long1, lat2, long2):

            # Convert latitude and longitude to
            # spherical coordinates in radians.
            degrees_to_radians = npy.pi/180.0

            # phi = 90 - latitude
            phi1 = (90.0 - lat1)*degrees_to_radians
            phi2 = (90.0 - lat2)*degrees_to_radians

            # theta = longitude
            theta1 = long1*degrees_to_radians
            theta2 = long2*degrees_to_radians

            # Compute spherical distance from spherical coordinates.

            # For two locations in spherical coordinates
            # (1, theta, phi) and (1, theta, phi)
            # cosine( arc length ) =
            #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
            # distance = rho * arc length
   
            cos = (npy.sin(phi1)*npy.sin(phi2)*npy.cos(theta1 - theta2) +
                   npy.cos(phi1)*npy.cos(phi2))
            arc = npy.arccos( cos )

            # Remember to multiply arc by the radius of the earth
            # in your favorite set of units to get length.
            arc[npy.where(self.maskrho == 0)] = 1.e36
            return arc

        def __call__(self,lonin,latin):
                a = self.distance_on_unit_sphere(latin, lonin, self.latrho, self.lonrho)

                #plt.figure() ; plt.pcolor(self.lonrho,self.latrho,a,vmin=0.,vmax=1.) ; plt.show()
                b = a.argmin()
                c = npy.unravel_index(b, self.lonrho.shape)
                print 'grid indexes are ', c , 'and in real life lon/lat = ', self.lonrho[c] , self.latrho[c]
                print 'and mask value is ', self.maskrho[c]
                return c

#------------------------------------------------------------------------
nwp = CCS_index('/t3/workdir/liz/MODELS/CTROMS/Inputs/coral_grd.nc')
#nwp = CCS_index('/t1/scratch/liz/Inputs/VIP-LD.HCo11T/Grid/VIP_grd_high_res_bathy_interp2.nc')

print 'Do not forget to add 1 to have fortran indexes'

print '#------------------------------------------------#'

lon_val = 120.0 + (58./60)
lat_val = 14.0 + (32./60)

print lon_val
print lat_val

fleet5 = nwp(lon_val,lat_val)
#fleet5 = nwp(120.966667,14.583333)

plt.figure()
plt.pcolor(nwp.lonvert,nwp.latvert,nwp.maskrho,cmap=cm.binary_r)
plt.plot(nwp.lonrho[fleet5],nwp.latrho[fleet5],'mo',linewidth=4)
plt.show()

