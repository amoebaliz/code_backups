import netCDF4 as nc
import numpy as npy
import matplotlib.pylab as plt
from matplotlib import cm as cm

class find_geo():

	def __init__(self,longrd,latgrd):
		self.lonrho = longrd
		self.latrho = latgrd
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
            #arc[npy.where(self.maskrho == 0)] = 1.e36
	    return arc

	def __call__(self,lonin,latin):
		a = self.distance_on_unit_sphere(latin, lonin, self.latrho, self.lonrho)

		#plt.figure() ; plt.pcolor(self.lonrho,self.latrho,a,vmin=0.,vmax=1.) ; plt.show()
		b = a.argmin()
		c = npy.unravel_index(b, self.lonrho.shape)
		print 'grid indexes are ', c , 'and in real life lon/lat = ', self.lonrho[c] , self.latrho[c]
		#print 'and mask value is ', self.maskrho[c]
		return c 

