#!/bin/bash

##################################################################
####              MODEL EVALUATION & ATTRIBUTION              ####
##################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTERNAL ARGUMENTS (PASSED FROM LOOP)

model_fnm=$1
gmst_fnm=$2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INTERNAL ARGUMENTS (COMMON TO ALL MODEL RUNS FOR SAME ANALYSIS)

id=62f4b5a82fde776a4c64f0ca33646aa0

distribution=gauss
fit_type=scale
include_event=on
return_period=10
ystart=1981
event_year=2022
restrain=0
gmst_past=-1.2
gmst_fut=0.8

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# store parameters as header line

header_file=tmp_headers.txt

info1=`echo "distribution fit_type include_event return_period obs_start event_year restrain gmst_past gmst_fut"`
info2=`echo "${distribution} ${fit_type} ${include_event} ${return_period} ${ystart} ${event_year} ${restrain} ${gmst_past} ${gmst_fut}"`

echo $info1 > $header_file
echo $info2 >> $header_file

####################################################################################################################################
####                                       SHOULDN'T NEED TO CHANGE ANYTHING AFTER THIS POINT                                   ####
####################################################################################################################################

# DERIVED VARIABLES

# get upload name from list of available files
model_line=`cat climexp_uploads.txt | grep $model_fnm`
upload_fnm=${model_line//*(/}
upload_fnm=${upload_fnm//)*/}

# get gmst name from list of available files
gmst_line=`cat climexp_uploads.txt | grep $gmst_fnm`
gmst_file=${gmst_line//*(/}
gmst_file=${gmst_file//)*/}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIXED VARIABLES (temporary files used to hold output)

eval_logfile=model_eval.log
attr_logfile=model_attr.log
proj_logfile=model_proj.log

####################################################################################################################################
####                                                 RUN ATTRIBUTION ON CLIMEXP                                                 ####
####################################################################################################################################

# RUN MODEL EVALUATION

data_script="https://climexp.knmi.nl/upload.cgi?STATION=$model_fnm&TYPE=$type&email=$id"
eval_script="https://climexp.knmi.nl/attribute.cgi?EMAIL=$id&STATION=$model_fnm&WMO=$upload_fnm&assume=$fit_type&begin=$ystart&biasrt=$return_period&ci=95&cov1=$gmst_past&dgt=80&end=$event_year&fit=$distribution&includelast=$include_event&key=right&restrain=$restrain&timeseries=./data/$gmst_file.1.$id.inf&type=attribute&year=$event_year"

curl ${data_script}
curl ${eval_script} > ${eval_logfile}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN MODEL ATTRIBUTION

attr_script="https://climexp.knmi.nl/attribute.cgi?EMAIL=$id&STATION=$model_fnm&WMO=$upload_fnm&assume=$fit_type&biasrt=$return_period&ci=95&cov1=$gmst_past&dgt=80&end=$event_year&fit=$distribution&includelast=$include_event&key=right&restrain=$restrain&timeseries=./data/$gmst_file.1.$id.inf&type=attribute&year=$event_year"

curl ${data_script}
curl ${attr_script} > ${attr_logfile}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN FUTURE ATTRIBUTION

# get value of current event
line=`cat $attr_logfile | grep "atr1"`; rv=${line//*(value/}; rv=`echo ${rv//)*/}`

proj_script="https://climexp.knmi.nl/attribute.cgi?EMAIL=$id&STATION=$model_fnm&WMO=$upload_fnm&assume=$fit_type&ci=95&cov1=$gmst_fut&dgt=80&end=2050&fit=$distribution&includelast=$include_event&key=right&restrain=$restrain&timeseries=./data/$gmst_file.1.$id.inf&type=attribute&xyear=$rv&year=$event_year"

curl ${data_script}
curl ${proj_script} > ${proj_logfile}


####################################################################################################################################
####                                                  EXTRACT RESULTS FROM LOG FILES                                            ####
####################################################################################################################################

## EVALUATION PARAMETERS - uncomment only the required lines, and modify headers accordingly

# # event_year
# line=`cat $eval_logfile | fgrep "&sigma;':" | fgrep "$event_year"`
# sigma_prime=`echo ${line#*</td><td>*${event_year}} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# dispersion
line=`cat $eval_logfile | fgrep "&sigma;/&mu;:"`
dispersion=`echo ${line#*</td><td>*} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# # shape parameter
# line=`cat $eval_logfile | grep "&xi;"`
# xi=`echo ${line#*</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# trend parameter
line=`cat $eval_logfile | grep "&alpha;"`
alpha=`echo ${line#*</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

eval_headers="disp_est disp_lower disp_upper alpha_est alpha_lower alpha_upper"
eval_params=`echo ${dispersion} ${alpha}`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ATTRIBUTION PARAMETERS

# return time & return value
line=`cat $attr_logfile | grep "atr1"`
rp=`echo ${line#*${event_year}</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# probability ratio
line=`cat $attr_logfile | grep "atra"`
pr=`echo ${line#*&nbsp;</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# change in intensity
line=`cat $attr_logfile | grep "change in intensity"`
DeltaI=`echo ${line#*diff*?</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

attr_headers="rv rp_est rp_lower rp_upper pr_est pr_lower pr_upper di_est di_lower di_upper"
attr_params=`echo ${rv} ${rp} ${pr} ${DeltaI}`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FUTURE ATTRIBUTION PARAMETERS
## Note that climexp computes attribution from future to present - currently this is handled in python postprocessing

# probability ratio
line=`cat $proj_logfile | grep "atra"`
pr_fut=`echo ${line#*&nbsp;</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# change in intensity
line=`cat $proj_logfile | grep "change in intensity"`
DeltaI_fut=`echo ${line#*diff*?</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# dispersion
line=`cat $proj_logfile | fgrep "&sigma;/&mu;:"`
dispersion_fut=`echo ${line#*</td><td>*} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# trend
line=`cat $proj_logfile | grep "&alpha;"`
alpha_fut=`echo ${line#*</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

proj_headers="f_pr_est f_pr_lower f_pr_upper f_di_est f_di_lower f_di_upper f_disp_est f_disp_lower f_disp_upper f_alpha_est f_alpha_lower f_alpha_upper"
proj_params=`echo ${pr_fut} ${DeltaI_fut} ${dispersion_fut} ${alpha_fut}`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write headers & parameters to file

header_line=`echo run gmst_file ${eval_headers} ${attr_headers} ${proj_headers}`
echo $header_line >> $header_file

model_line=`echo ${model_fnm} ${gmst_fnm} ${eval_params} ${attr_params} ${proj_params}`
echo $model_line >> tmp_results.txt

