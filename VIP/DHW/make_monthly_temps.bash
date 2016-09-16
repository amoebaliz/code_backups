#!/bin/bash

DIR=/t3/workdir/liz/MODELS/VIP/Runs/VIP-LD.HCo13T/
dest_dir=/t3/workdir/liz/MODELS/VIP/PostProc/monthly_temps/ 
filstr=VIP-LD.HCo13T_avg_

# clean up
rm ${dest_dir}tmp_*

# loop over years and months
for yr in {1996..1999} ; 
do

    echo $yr
    for mon in {1..12} ;
    do
        
        echo $mon
        monstr=$(printf %02d $mon)
        ncra $DIR$yr/$filstr$yr-$monstr-* ${dest_dir}tmp_$monstr.nc
    done
    ncrcat ${dest_dir}tmp* ${dest_dir}monthly_temps_$yr.nc
    rm ${dest_dir}tmp* 

done

