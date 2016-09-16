#!/bin/bash


# loop over all files
for num in {17901..18262} ; 
#for num in {17900..17900} ;
do
 filename=coral_avg_${num}_VIP.nc
 vel_file=./ctroms_vel/$filename
 temp_file=./ctroms_temp/$filename
 salt_file=./ctroms_salt/$filename

 # append temp to vel
 ncks -A -v temp $temp_file $vel_file

 # append salt to vel
 ncks -A -v salt $salt_file $vel_file

done

