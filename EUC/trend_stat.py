def trend_stat(M):
	from scipy import stats
	import numpy as np
# INPUT: 
# M is the input time series <nt x 1>
# siglev is the significance level (e.g., 95)

# OUTPUT:
# trend is the linear trend
# plusminus is the value for trend #/- plusminus
# sig is 1 for significant, 0 for not significant
# trend_pom is trend expressed as precent of mean (pom)

# NOTES:

# Use MATLAB's linear function to calculate a trend.

# Does not explicitely take into account degrees of freedom, so significance may not be valid if time series is smoothed/filtered. 

# Should still work if there are NaNs in the input time series. 

# The units of trend are the units of the value in the input time series per unit increment of time series

# Time spacing must be linear.
# UPDATE: EJD - added if statememt to deal w/ all NaN timeseries


	X = np.ones((len(M),2))
	X[:,1]=np.arange(1,len(M)+1)
	
	if all(np.isnan(M)):
		trend = NaN
		sig = 0
		plusminus = 0
		trend_pom = 0
	else:
	
		yM = M
		print M.shape
		print X.shape
		slope, intercept, r_value, p_value, std_err = stats.linregress(yM,X)
		trend = slope

		sig = 0
		plusminus = std_err
	
	if abs(trend) > abs(plusminus):
		sig = 1

	 
	trend_pom = 100*trend/mean(M)  
	return trend_pom
	return sig




