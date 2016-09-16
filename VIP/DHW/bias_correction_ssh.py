#!/usr/bin/env python

import numpy as np
from lib_bias_correction_ssh import bias_correction_ssh as bc

dir_differences = '/Volumes/P4/workdir/raphael/analysis_CCS1-Cobalt/bias_correction/1_compare_woa13/differences/'
mask_file = './mask_GFDL.nc'

dir_mom_historical_physics = '/Volumes/P4/workdir/raphael/RUN_GFDL/processed/ocean_physics/daily/'
dirout_mom_historical_physics = '/Volumes/P4/workdir/raphael/analysis_CCS1-Cobalt/bias_correction/2_corrected_historical/physics/daily/'
mom_runname = '_GFDL-ESM2Mnew_'
first_year = 1951
last_year  = 2005

#--------------------------------------------------------------------------------------

days_in_month = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

for year in np.arange(first_year,last_year+1):
	print 'Working on year ', year
	yyyy=str(year).zfill(4)
	nn_day = 0
	for month in np.arange(12):
		mm=str(month+1).zfill(2)
		for day in np.arange(days_in_month[month]):
			dd=str(day+1).zfill(2)
			datetag='y' + yyyy + 'm' + mm + 'd' + dd
			nn_day = nn_day + 1
			var = 'eta_t'
			file_mom_in = 'eta' + mom_runname + datetag + '.nc'
			diff_file = dir_differences + 'monthly_difference_' + var + '_ESM2M-SODA.nc'

			bias = bc(diff_file,mask_file,dir_mom_historical_physics,dirout_mom_historical_physics)
			bias(var,file_mom_in,nn_day)

