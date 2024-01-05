#!/bin/bash

# script to subset each subfile of the large ensemble (concatenation done separately)
module load cdo

varnm=snd
for gcm in `seq 1 5`; do 
    for rcm in `seq 1 7`; do 
        
        fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/00_data-for-fwi/canesm-canrcm/subfiles/${varnm}_NAM-44_CCCma-CanESM2_historical-r${gcm}_r${rcm}i1p1_CCCma-CanRCM4_r2_*.nc`;
        new_fnm=canesm-canrcm/${varnm}_NAM-44_CCCma-CanESM2_historical-r${gcm}_r${rcm}i1p1_CCCma-CanRCM4_r2.nc;
        cdo cat $fl $new_fnm;
        
    done;
done