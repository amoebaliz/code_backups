import pyroms
import numpy as np

class cfl():

	def __init__(self,domain):
		# get my roms grid
		self.grd = pyroms.grid.get_ROMS_grid(domain)
		# earth gravity in m.s-2
		self.g = 9.81
		return None

	def __call__(self,dt,ndtfast):
		# store input args
		self.my_dt = float(dt)
		self.my_ndtfast = float(ndtfast)
		# compute gravity waves speed
		self.gravity_waves_max_speed()
		# find minimum distance in domain
		self.min_distance_grid()
		# compute speed of gravity waves
		self.max_barotropic_dt()
		# check if values are ok
		self.check_my_values()
		return None

	def gravity_waves_max_speed(self):
		# maximum depth in domain
		hmax = self.grd.vgrid.h.max()
		# speed of gravity waves at max depth = sqrt( g * hmax)
		self.gw_speed = np.sqrt( self.g * hmax )
		print '>>> Max Speed of gravity waves is', self.gw_speed, 'm.s-1'
		return None

	def min_distance_grid(self):
		# read x,y arrays from grid
		x_vert = self.grd.hgrid.x_vert.copy()
		y_vert = self.grd.hgrid.y_vert.copy()
		# compute distances in meters
		dx = np.abs(x_vert[:,1:] - x_vert[:,:-1])
		dy = np.abs(y_vert[1:,:] - y_vert[:-1,:])
		# take minimum for both directions
		dxmin = dx.min()
		dymin = dy.min()
		# and then the minimum of those two
		self.dxymin = np.minimum(dxmin,dymin)
		print '>>> Minimum distance between points is', self.dxymin, 'm'
		return None
		
	def max_barotropic_dt(self):
		# my theoritical barotropic timestep dt < dx / sqrt(g.h)
		self.max_bdt = self.dxymin / self.gw_speed
		print '>>> Theoritical maximum barotropic timestep is', self.max_bdt , 's'
		return None

	def check_my_values(self):
		if self.my_ndtfast < 15:
			print 'value too low for ndtfast'
		if self.my_ndtfast > 55:
			print 'value too high for ndtfast'

		# compute my barotrpoic timestep
		my_bdt = self.my_dt / self.my_ndtfast
		# check against theoritical value
		if my_bdt <= self.max_bdt:
			print '>>> dt and ndtfast seem ok'
			print '>>> dt / ndtfast =', my_bdt , '<=', self.max_bdt
		else:
			print '>>> dt and ndtfast NOT ok'
			print '>>> dt / ndtfast =', my_bdt , '>', self.max_bdt
		return None
		

#--------------- example ------------------

if __name__ == "__main__":

	import libcfl
	mycfl = libcfl.cfl('CCS2')
	mycfl(200,30)
