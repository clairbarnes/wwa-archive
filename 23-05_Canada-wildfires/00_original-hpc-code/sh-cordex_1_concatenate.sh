#!/bin/bash

# script to concatenate subfiles of CORDEX runs
module load cdo

for varnm in tas sfcWind snw hurs; do
    for mdl in HadGEM2 MPI-ESM NorESM; do 
        
        fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/cordex/${varnm}_*${mdl}*.nc`;
        
        fnm_root=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/cordex/${varnm}_*${mdl}*.nc | head -1`;
        fnm_root=${fnm_root##*/};
        new_fnm=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/00_data-for-fwi/cordex/${fnm_root/_1970*/.nc};
        cdo cat $fl $new_fnm; 
    done
done