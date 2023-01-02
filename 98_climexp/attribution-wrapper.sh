
####################################################################################################################################
####                                  MODEL EVALUATION & ATTRIBUTION USING THE CLIMATE EXPLORER                                 ####
####################################################################################################################################

# filename for final results
results_file=results_lower-niger.txt

# list of model:gmst pairs to loop over, separated with spaces
model_list=p_lower-niger_7daymax_NorESM1-M_r1_RegCM4-7:tas_Amon_NorESM1-M_rcp85_r1i1p1
#`echo p_lower-niger_7daymax_NorESM1-M_r1_RegCM4-7:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_EC-EARTH_r12_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_HadGEM2-ES_r1_RegCM4-3:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_NorESM1-M_r1_CCLM5-0-15:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MPI-ESM-LR_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MPI-ESM-MR_r1_RegCM4-7:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MPI-ESM-MR_r1_RegCM4-3:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_IPSL-CM5A-LR_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_HadGEM2-ES_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_HadGEM2-ES_r1_REMO2015:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lower-niger_7daymax_MIROC5_r1_REMO2009:tas_Amon_NorESM1-M_rcp85_r1i1p1`

# list of parameter=value pairs to use in study, separated with spaces
parameter_list=`echo distribution=gev fit_type=scale return_period=5 obs_start=1981 event_year=2022 lower_tail= restrain=0.4 include_event=on gmst_past=-1.2 gmst_fut=0.8 confint=95`

# id of user that uploaded the data
user_id=62f4b5a82fde776a4c64f0ca33646aa0


##############################################################################################################################
####                                  SHOULDN'T NEED TO CHANGE ANYTHING AFTER THIS LINE                                   ####
##############################################################################################################################


## DOWNLOAD LIST OF AVAILABLE DATA ON CLIMATE EXPLORER

curl https://climexp.knmi.nl/userseries.cgi?id=$user_id > climexp_uploads.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## CALL CLIMATE EXPLORER TO RUN ATTRIBUTION FOR EACH MODEL IN TURN

# loop over the models and run the attribution
for m in $model_list
do
    # split model string into model data & GMST string
    model_fnm=`echo ${m/:*/}`
    gmst_fnm=`echo ${m/*:/}`
    
    bash attribution-loop.sh $model_fnm $gmst_fnm "${parameter_list}" $user_id
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## CONCATENATE HEADERS & RESULTS INTO A SINGLE FILE

cat tmp_headers.txt > $results_file
cat tmp_results.txt >> $results_file

# remove temporary files
rm tmp_headers.txt tmp_results.txt model_eval.log model_attr.log model_proj.log
