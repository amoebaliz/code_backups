import numpy as np
import netCDF4 as nc 
import matplotlib.pylab as plt
import static_instab_correction as static

class bias_correction():

	def __init__(self,diff_file,maskfile,dirin,dirout):
		# path to obs files,...	
		self.days_in_month = np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
		self.diff_file = diff_file
		self.dirin  = dirin
		self.dirout = dirout
		self.spval = -1.0e+20
		self.maskfile = maskfile
		return None

	def __call__(self,varname,filename,nn_day):
		self.varname = varname
		self.filename = self.dirin + filename
		self.nn_day = nn_day

		self.fileout = self.dirout + filename.replace('ESM2Mnew','ESM2MnewBiasCorr')

		my_weights = self.time_interp_weights(self.nn_day)
		my_bias = self.interp_bias_to_day(my_weights)


		mynewdata = self.correct_the_bias(my_bias)

		self.clone_ncfile(self.filename,self.fileout)
		self.correct_ncfile(self.fileout,self.varname,mynewdata)

		return None

	def correct_the_bias(self,my_bias):
		original_field = self.readnc(self.filename,self.varname)
		depth = self.readnc(self.filename,'st_ocean')
		nz,ny,nx = original_field.shape

		mask = self.readnc(self.maskfile,'lsm')	

		# level where there is no observations 
		if self.varname in ['temp','salt']:
			self.depth_of_no_data = 2000.
			self.first_level_nodata = 39 # 40 in fortran
			self.last_level_data = 39 -1 # 39 in fortran
			self.level_check_static = 34 # 1000 meters

		# in the upper part, we remove the mean bias
		upper_corrected = original_field - my_bias
		upper_corrected[self.first_level_nodata:,:,:] = 0.

		# the lower part is more complicated
		# we add the original shape to last corrected value
		lower_corrected = np.zeros((nz,ny,nx))
		delta_bottom    = np.zeros((nz,ny,nx))

		original_bottom = original_field.copy()

		delta_bottom = original_bottom[:,:,:] - original_field[self.last_level_data,:,:]
		lower_corrected = upper_corrected[self.last_level_data,:,:] + delta_bottom[:,:,:]
		
		lower_corrected[:self.first_level_nodata,:,:] = 0.

		# add correction for the bottom
		all_corrected = upper_corrected + lower_corrected

		# correct for static instab
		if self.varname == 'temp':
			all_corrected = static.correct_static_temperature(all_corrected,self.level_check_static)
		elif self.varname == 'salt':
			all_corrected = static.correct_static_salinity(all_corrected,self.level_check_static)
		else:
			pass

		# correct for stupid values
		if self.varname == 'temp':
			# limit min temperature to -1.8 C
			all_corrected[np.where(all_corrected < -1.8 + 273.15)] = -1.8 + 273.15
		elif self.varname == 'salt':
			all_corrected[np.where(all_corrected < 0.)] = 0.
		else:
			pass
	
		# masked values with original mask
		all_corrected[np.where(mask == 0)] = self.spval
		return all_corrected

	def interp_bias_to_day(self,my_weights):
		# time interpolation of bias file
		bias_monthly = self.readnc(self.diff_file,self.varname)
		nt,nz,ny,nx = bias_monthly.shape
		if ( nt != 12):
			print 'Error : number of frames' 

		bias_thisday = np.zeros((nz,ny,nx))
		for kt in np.arange(12):
			bias_thisday = bias_thisday + bias_monthly[kt,:,:,:] * my_weights[kt]

		return bias_thisday

	def time_interp_weights(self,nday):
		''' compute weights to apply to a monthly data file
		nday is numbers of days since begining of the year '''
		# init to zero
		weights = np.zeros((12))
		# cumul of days at the middle of the month
		cumul_days = self.days_in_month.cumsum() - (self.days_in_month / 2.)
		# find which month is the lowest bound of the interval
		# uses the property of python that -1 is last index
		nmonth = -1
		for km in np.arange(12):
			if ( nday >= cumul_days[km] ):
				nmonth = nmonth + 1
		# distance from current day to lower bound of interval
		dist_to_lower_bound = np.mod(nday - cumul_days[nmonth],365)
		# distance between lower and upper bound (uses periodicity properties for nmonths and ndays)
		interval = np.mod(cumul_days[np.mod(nmonth+1,12)] - cumul_days[np.mod(nmonth,12)],365)
		# weight for lower and upper bound
		wgt_1 = 1.0 - (float(dist_to_lower_bound) / interval )
		wgt_2 = 1.0 - wgt_1
		# change only values for lower and upper bound (all other remain zero)
		weights[np.mod(nmonth,12)] = wgt_1
		weights[np.mod(nmonth+1,12)] = wgt_2
		return weights

	def verif_time_interp_algo(self):
		''' test for the weight algo '''
		weights_all = np.empty((365,12))
		for kday in np.arange(365):
			weights_all[kday,:] = self.time_interp_weights(kday)

		plt.figure()
		plt.plot(weights_all[:,:])
		plt.xticks(self.days_in_month.cumsum(),self.days_in_month.cumsum())
		plt.show()
		return None

        def readnc(self,filein,varin):
                ''' read data from netcdf '''
		fid = nc.Dataset(filein,'r')
                out = fid.variables[varin][:].squeeze()
		fid.close()
                return out

        def clone_ncfile(self,filein,fileout):
		fidin  = nc.Dataset(filein,'r')
                fidout = nc.Dataset(fileout,'w',format='NETCDF3_CLASSIC')

		fidout.description = 'ESM2M bias corrected (raphael.dussin@gmail.com)'
		#Copy dimensions
		for dname, the_dim in fidin.dimensions.iteritems():
			fidout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

		# Copy variables
		for v_name, varin in fidin.variables.iteritems():
			outVar = fidout.createVariable(v_name, varin.datatype, varin.dimensions)
			# Copy variable attributes
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
			outVar[:] = varin[:]

		fidin.close()
		fidout.close()
                return None

	def correct_ncfile(self,ncfile,varname,data):
		fid = nc.Dataset(ncfile,'a')
		ovariable = fid.variables[varname]
		ovariable[0,:,:,:] = data
		fid.close()
		return None
