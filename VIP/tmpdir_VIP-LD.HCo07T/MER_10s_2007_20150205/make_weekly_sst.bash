#!/bin/bash

day=7
for nweek in {1..50} ; 
do
 printf -v DAY "%04d" $day
 printf -v NWEEK "%02d" $nweek
 echo $DAY 
 ncra -O -n 7,4,1 -d s_rho,-1 -v temp VIP-LD.HCo07T_avg_$DAY.nc tmp_week_$NWEEK.nc
 day=$((day+7))
done
ncrcat -O tmp_week_*.nc weekly_vip_sst.nc
