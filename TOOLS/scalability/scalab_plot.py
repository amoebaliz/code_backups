import numpy as npy
import os
import matplotlib.pylab as plt

class scalability_test():

	def __init__(self,listlogs,indict):
		self.listlogs = listlogs
		self.NN = len(self.listlogs)
		# keys from input dict
		self.machine    = indict['machine']
		self.experiment = indict['experiment']
		self.compiler   = indict['compiler']
		self.timesteps  = indict['timesteps']
		self.nbpoints   = indict['nbpoints']
		# init sizes from number of logs
		self.ncores   = npy.zeros((self.NN))
		self.wallmax  = npy.zeros((self.NN))
		self.wallmin  = npy.zeros((self.NN))
		self.wallmean = npy.zeros((self.NN))
		# read logs
		self.read_logs()
		# compute speedup/efficiency
		self.compute_speedup()
		self.compute_efficiency()
		# plots
		self.plot_ellapsed()
		self.plot_speedup()
		#self.plot_efficiency()
		return None

	def read_logs(self):
		for nl in npy.arange(self.NN):
			self.ncores[nl] = self.get_core_info(self.listlogs[nl])
			self.wallmin[nl], self.wallmean[nl], self.wallmax[nl] = self.get_walltime(self.listlogs[nl])
		return None

	def compute_speedup(self):
		indref = self.ncores.argmin()
		self.ncmin = self.ncores.min()
		self.speedup = self.wallmax[indref] / self.wallmax
		return None
	
	def compute_efficiency(self):
		indref = self.ncores.argmin()
		self.efficient = self.wallmax[indref] / ( self.ncores * self.wallmax )
		return None

	def get_core_info(self,roms_logfile):
		#''' Get number of points, number of cores and layout
		#    from ROMS log file '''
		command = 'grep Tiling ' + roms_logfile
		out = os.popen(command).readlines()
		#npoints = out[0].split()[3].replace(',','')
		ncores  = out[0].split()[6].replace(',','')
		ncores  = float(ncores)
		#tiling  = out[0].split()[8]
		#return npoints, ncores, tiling
		return ncores

	def get_walltime(self,roms_logfile):
		#''' Get the walltime from ROMS log file '''
		command = "grep 'Node   #' " + roms_logfile + " | awk '{ print $NF }' "
		out = os.popen(command).readlines()
		wall = []
		for wt in out:
			wall.append(float(wt.replace('\n','')))
		wall     = npy.array(wall)
		wallmin  = wall.min()
		wallmean = wall.mean()
		wallmax  = wall.max()
		return wallmin, wallmean, wallmax

	def plot_ellapsed(self):
		plt.figure(figsize=[8.,8.])
		plt.plot(self.ncores,self.wallmax/60.,'k')
		plt.grid()
		plt.xlabel('Number of Cores')
		plt.ylabel('Elapsed time (min)')
		mytitle = self.timesteps + ' timesteps of ' + self.experiment + ' run on machine ' + self.machine + \
		'\nNumber of points : ' + self.nbpoints + \
	        '\nROMS compiled with ' + self.compiler
		plt.title(mytitle)
		myfilename = 'Ellapsed_' + self.experiment + '_' + self.machine + '.' + self.compiler + \
		'_' + self.timesteps + 'steps' + '.png'
		plt.savefig(myfilename)
		return None

	def plot_speedup(self):
		plt.figure()
		plt.plot(self.ncores,self.ncores/self.ncmin,'k--',self.ncores,self.speedup,'k')
		plt.grid()
		plt.xlabel('Number of Cores')
		plt.ylabel('Speedup')
		mytitle = self.timesteps + ' timesteps of ' + self.experiment + ' run on machine ' + self.machine + \
		'\nNumber of points : ' + self.nbpoints + \
	        '\nROMS compiled with ' + self.compiler
		plt.title(mytitle)
		myfilename = 'Speedup_' + self.experiment + '_' + self.machine + '.' + self.compiler + \
		'_' + self.timesteps + 'steps' + '.png'
		plt.savefig(myfilename)
		return None

	def plot_efficiency(self):
		plt.figure()
		plt.plot(self.ncores,self.efficient,'k')
		plt.grid()
		plt.xlabel('Number of Cores')
		plt.ylabel('Efficiency')
		plt.title('Triton : ')
		plt.savefig('./efficiency.png')
		return None




