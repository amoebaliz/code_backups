import datetime as dt


#creating the relativedelta to skip months 
def addmonth(t_ref,n):
	
	# decompose t_ref 
	year0 = t_ref[0]
	month0 = t_ref[1]
	day0 = t_ref[2]
	
	# first figure out how many years	
	year1 = year0 + n//12
	month1 = month0 + n%12

		if month1 > 12:
			year1 = year1 + 1
			month1 = month1 - 12	
	
	# deal with days
	
	deci_month = month1 % 1
	if deci_month == 0:
		day1=day0
	else:
		# look up the month and choose a day based on the fraction
		day1 = np.ceil(deci_month * cal.itermonthdays(year1,month1))

	retur dt.date(year1,month1,day1)


