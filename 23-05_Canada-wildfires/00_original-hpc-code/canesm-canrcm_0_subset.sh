#!/bin/bash

# script to subset each subfile of the large ensemble (concatenation done separately)
module load cdo

# temperature (only 18:00)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/can_lens/*/tas_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/00_data-for-fwi/canesm-canrcm/subfiles/${file_in##*/};
    cdo -s selhour,18 -sellonlatbox,280,297,47,59 $file_in $file_out;
done
echo "tas complete"

# relative humidity (only 18:00)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/can_lens/*/hurs_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/00_data-for-fwi/canesm-canrcm/subfiles/${file_in##*/};
    cdo -s selhour,18 -sellonlatbox,280,297,47,59 $file_in $file_out;
done
echo "hurs complete"

# wind speed (only 18:00)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/can_lens/*/sfcWind_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/00_data-for-fwi/canesm-canrcm/subfiles/${file_in##*/};
    cdo -s selhour,18 -sellonlatbox,280,297,47,59 $file_in $file_out;
done
echo "sfcWind complete"

# precip (keep all timesteps)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/can_lens/*/pr_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/00_data-for-fwi/canesm-canrcm/subfiles/${file_in##*/};
    cdo -s sellonlatbox,280,297,47,59 $file_in $file_out;
done
echo "pr complete"

# snow depth (keep all timesteps)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/can_lens/*/snd_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/00_data-for-fwi/canesm-canrcm/subfiles/${file_in##*/};
    cdo -s sellonlatbox,280,297,47,59 $file_in $file_out;
done;
echo "snd complete"
