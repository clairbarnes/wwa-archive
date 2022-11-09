
model_list=`echo p_lower-niger_7daymax_NorESM1-M_r1_RegCM4-7:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_EC-EARTH_r12_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_HadGEM2-ES_r1_RegCM4-3:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_NorESM1-M_r1_CCLM5-0-15:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MPI-ESM-LR_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MPI-ESM-MR_r1_RegCM4-7:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MPI-ESM-MR_r1_RegCM4-3:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_IPSL-CM5A-LR_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_HadGEM2-ES_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_HadGEM2-ES_r1_REMO2015:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MIROC5_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1`

results_file=results_lower-niger.txt
attr_script=climexp_attribution_lower-niger.sh


##############################################################################################################################
####                                  SHOULDN'T NEED TO CHANGE ANYTHING AFTER THIS LINE                                   ####
##############################################################################################################################


## DOWNLOAD LIST OF AVAILABLE DATA ON CLIMATE EXPLORER

curl https://climexp.knmi.nl/userseries.cgi?id=62f4b5a82fde776a4c64f0ca33646aa0 > climexp_uploads.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## CALL CLIMATE EXPLORER TO RUN ATTRIBUTION FOR EACH MODEL IN TURN

# loop over the models and run the attribution
for m in $model_list
do
    # split model string into model data & GMST string
    model_fnm=`echo ${m/:*/}`
    gmst_fnm=`echo ${m/*:/}`
    
    bash $attr_script $model_fnm $gmst_fnm
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## CONCATENATE HEADERS & RESULTS INTO A SINGLE FILE

cat tmp_headers.txt > $results_file
cat tmp_results.txt >> $results_file

# remove temporary files
rm tmp_headers.txt tmp_results.txt model_eval.log model_attr.log model_proj.log