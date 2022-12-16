
####################################################################################################################################
####                                  MODEL EVALUATION & ATTRIBUTION USING THE CLIMATE EXPLORER                                 ####
####################################################################################################################################

# filename for final results
results_file=cx_station-fit.txt

# list of model:gmst pairs to loop over, separated with spaces
model_list=`echo t_AR-PY_tasmax-7day_station_87129 t_AR-PY_tasmax-7day_station_87078 t_AR-PY_tasmax-7day_station_87065 t_AR-PY_tasmax-7day_station_87585 t_AR-PY_tasmax-7day_station_87345 t_AR-PY_tasmax-7day_station_87244 t_AR-PY_tasmax-7day_station_87623 t_AR-PY_tasmax-7day_station_87480`

# list of parameter=value pairs to use in study, separated with spaces
parameter_list=`echo distribution=gev fit_type=shift event_year=2022 lower_tail= restrain=0.4 include_event=on gmst_past=-1.2 gmst_fut=0.8 confint=95`

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
    
    bash check-obs-fit_loop.sh $model_fnm $gmst_fnm "${parameter_list}" $user_id
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## CONCATENATE HEADERS & RESULTS INTO A SINGLE FILE

cat tmp_headers.txt > $results_file
cat tmp_results.txt >> $results_file

# remove temporary files
rm tmp_headers.txt tmp_results.txt model_eval.log climexp_uploads.txt
