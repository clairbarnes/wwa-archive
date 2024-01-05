#!/bin/bash

# script to subset each subfile of the large ensemble (concatenation done separately)
# Chaining fails for some reason, so using an intermediate temporary file
module load cdo

tmp_file=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/tmp.nc

# # temperature (only 16:30)
# fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/synda/data/*-*/*/*/tas/tas_*.nc`
# for file_in in $fl; do
#     file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/cordex/${file_in##*/};
#     cdo -s selindexbox,197,240,135,195 $file_in $tmp_file;
#     cdo -s selhour,16 $tmp_file $file_out;
# done
# echo "tas complete"

# relative humidity (only 16:30)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/synda/data/*-*/*/*/hurs/hurs_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/cordex/${file_in##*/};
    cdo -s selindexbox,197,240,135,195 $file_in $tmp_file;
    cdo -s selhour,16 $tmp_file $file_out;
done
echo "hurs complete"

# wind speed (only 16:30)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/synda/data/*-*/*/*/sfcWind/sfcWind_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/cordex/${file_in##*/};
    cdo -s selindexbox,197,240,135,195 $file_in $tmp_file;
    cdo -s selhour,16 $tmp_file $file_out;
done
echo "sfcWind complete"

# precip (all)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/synda/data/*-*/*/*/pr/pr_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/cordex/${file_in##*/};
    cdo -s selindexbox,197,240,135,195 $file_in $file_out;
done
echo "pr complete"

# snow (all)
fl=`ls /rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/synda/data/*-*/*/*/snw/snw_*.nc`
for file_in in $fl; do
    file_out=/rds/general/user/cb2714/home/00_WWA_project_folder/ephemeral/canada_fwi/99_processing/cordex/${file_in##*/};
    cdo -s selindexbox,197,240,135,195 $file_in $file_out;
done;
echo "snw complete"