
####################################################################################################################################
####                                  MODEL EVALUATION & ATTRIBUTION USING THE CLIMATE EXPLORER                                 ####
####################################################################################################################################

# filename for final results
results_file=cx-results_highresmip.txt

# list of model:gmst pairs to loop over, separated with spaces
model_list=p_rx2d_NZ_CMCC-CM2-VHR4 #`echo p_rx2d_NZ_MPI-ESM1-2-XR p_rx2d_NZ_EC-Earth3P-HR p_rx2d_NZ_CNRM-CM6-1-HR p_rx2d_NZ_HadGEM3-GC31-HM  p_rx2d_NZ_HadGEM3-GC31-MM`

# list of parameter=value pairs to use in study, separated with spaces
parameter_list=`echo distribution=gev fit_type=scale return_period=10 obs_start=1979 event_year=2023 lower_tail= restrain=0.4 include_event=on gmst_past=-1.2 gmst_fut=0.8 confint=95`

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
    
    bash cx-attribution-loop-HighResMIP.sh $model_fnm "${parameter_list}" $user_id
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## CONCATENATE HEADERS & RESULTS INTO A SINGLE FILE

cat tmp_headers.txt > $results_file
cat tmp_results.txt >> $results_file

# remove temporary files
rm tmp_headers.txt tmp_results.txt model_eval.log model_attr.log model_proj.log climexp_uploads.txt
