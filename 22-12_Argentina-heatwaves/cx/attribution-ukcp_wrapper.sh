
####################################################################################################################################
####                                  MODEL EVALUATION & ATTRIBUTION USING THE CLIMATE EXPLORER                                 ####
####################################################################################################################################

# filename for final results
results_file=cx_attribution-ukcp.txt

# list of model:gmst pairs to loop over, separated with spaces
model_list=`echo t_AR-PY_tasmax-7day_UKCP-land-gcm_01:t_ukcp-land-gcm_gsat-aw_sm_01 t_AR-PY_tasmax-7day_UKCP-land-gcm_02:t_ukcp-land-gcm_gsat-aw_sm_02 t_AR-PY_tasmax-7day_UKCP-land-gcm_03:t_ukcp-land-gcm_gsat-aw_sm_03 t_AR-PY_tasmax-7day_UKCP-land-gcm_04:t_ukcp-land-gcm_gsat-aw_sm_04 t_AR-PY_tasmax-7day_UKCP-land-gcm_05:t_ukcp-land-gcm_gsat-aw_sm_05 t_AR-PY_tasmax-7day_UKCP-land-gcm_06:t_ukcp-land-gcm_gsat-aw_sm_06 t_AR-PY_tasmax-7day_UKCP-land-gcm_07:t_ukcp-land-gcm_gsat-aw_sm_07 t_AR-PY_tasmax-7day_UKCP-land-gcm_08:t_ukcp-land-gcm_gsat-aw_sm_08 t_AR-PY_tasmax-7day_UKCP-land-gcm_09:t_ukcp-land-gcm_gsat-aw_sm_09 t_AR-PY_tasmax-7day_UKCP-land-gcm_10:t_ukcp-land-gcm_gsat-aw_sm_10 t_AR-PY_tasmax-7day_UKCP-land-gcm_11:t_ukcp-land-gcm_gsat-aw_sm_11 t_AR-PY_tasmax-7day_UKCP-land-gcm_12:t_ukcp-land-gcm_gsat-aw_sm_12 t_AR-PY_tasmax-7day_UKCP-land-gcm_13:t_ukcp-land-gcm_gsat-aw_sm_13 t_AR-PY_tasmax-7day_UKCP-land-gcm_14:t_ukcp-land-gcm_gsat-aw_sm_14 t_AR-PY_tasmax-7day_UKCP-land-gcm_15:t_ukcp-land-gcm_gsat-aw_sm_15`

# list of parameter=value pairs to use in study, separated with spaces
parameter_list=`echo distribution=gev fit_type=shift return_period=20 obs_start=1950 event_year=2022 lower_tail= restrain=0.4 include_event=on gmst_past=-1.2 gmst_fut=0.8 confint=95`

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
    
    bash attribution-ukcp_loop.sh $model_fnm $gmst_fnm "${parameter_list}" $user_id
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## CONCATENATE HEADERS & RESULTS INTO A SINGLE FILE

cat tmp_headers.txt > $results_file
cat tmp_results.txt >> $results_file

# remove temporary files
rm tmp_headers.txt tmp_results.txt model_eval.log model_attr.log model_proj.log climexp_uploads.txt
