#!/bin/bash

####################################################################################################################################
####                                  MODEL EVALUATION & ATTRIBUTION USING THE CLIMATE EXPLORER                                 ####
####################################################################################################################################

# STILL TO DO:

    # are extra parameters needed for other distributions?
    
####################################################################################################################################
# EXTERNAL ARGUMENTS (PASSED FROM WRAPPER)

model_fnm=$1
gmst_fnm=$2
par_list=$3
id=$4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# DERIVED VARIABLES

# split list of parameters into individual variables
for par in $par_list; do eval ${par}; done

# get upload name from list of available files
model_line=`cat climexp_uploads.txt | grep kalelink | grep $model_fnm`
upload_fnm=${model_line//*(/}
upload_fnm=${upload_fnm//)*/}

# get gmst name from list of available files
gmst_file=gmst

# don't change sign unless 'lower_tail' is passed as parameter
chsign="&changesign=on"
if [[ -z "$lower_tail" ]]; then chsign=""; fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FIXED VARIABLES (temporary files used to hold output, will be deleted after processing)

eval_logfile=model_eval.log
attr_logfile=model_attr.log
proj_logfile=model_proj.log
tmp_header_file=tmp_headers.txt
tmp_result_file=tmp_results.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# STORE PARAMETERS AS DUMMY HEADER LINE

echo "# ${par_list}" > $tmp_header_file


####################################################################################################################################
####                                                 RUN ATTRIBUTION ON CLIMEXP                                                 ####
####################################################################################################################################

# RUN MODEL EVALUATION

data_script="https://climexp.knmi.nl/upload.cgi?STATION=$model_fnm&email=$id"
eval_script="https://climexp.knmi.nl/attribute.cgi?EMAIL=$id&STATION=$model_fnm&WMO=$upload_fnm&assume=$fit_type&biasrt=$return_period&ci=$confint&cov1=$gmst_past&dgt=80&fit=$distribution&includelast=$include_event&restrain=$restrain&timeseries=$gmst_file&type=attribute&year=$event_year"$chsign  
  
curl ${data_script}
curl ${eval_script} > ${eval_logfile}

# get value of current event (needed for future attribution)
line=`cat $eval_logfile | grep "atr1"`; rv=${line//*(value/}
rv=`echo ${rv//)*/}`


####################################################################################################################################
####                                                  EXTRACT RESULTS FROM LOG FILES                                            ####
####################################################################################################################################

## EVALUATION PARAMETERS - uncomment only the required lines, and modify headers accordingly

# sigma_prime
line=`cat $eval_logfile | fgrep "&sigma;':" | fgrep "$event_year"`
sigma_prime=`echo ${line#*</td><td>*${event_year}} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# dispersion
line=`cat $eval_logfile | fgrep "&sigma;/&mu;:"`
dispersion=`echo ${line#*</td><td>*} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`
disp_header="disp_est disp_lower disp_upper"; if [[ -z "$dispersion" ]]; then disp_header=""; fi

# shape parameter (remove header if parameter is empty)
line=`cat $eval_logfile | grep "&xi;"`
xi=`echo ${line#*</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`
xi_header="xi_est xi_lower xi_upper"; if [[ -z "$xi" ]]; then xi_header=""; fi

# trend parameter
line=`cat $eval_logfile | grep "&alpha;"`
alpha=`echo ${line#*</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

eval_headers=`echo "sigma_est sigma_lower sigma_upper ${disp_header} ${xi_header} alpha_est alpha_lower alpha_upper"`
eval_params=`echo ${sigma_prime} ${dispersion} ${xi} ${alpha}`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## ATTRIBUTION PARAMETERS


# return time & return value
line=`cat $eval_logfile | grep "atr1"`
rp=`echo ${line#*${event_year}</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# probability ratio
line=`cat $eval_logfile | grep "atra"`
pr=`echo ${line#*&nbsp;</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

# change in intensity
line=`cat $eval_logfile | grep "change in intensity"`
DeltaI=`echo ${line#*diff*</td><td>} | sed -e 's_</td><td>__g' | sed -e 's_</td></tr>__g' | sed -e 's_\.\.\.__g' | sed -e 's/&infin;/1e6/g'`

attr_headers="rv rp_est rp_lower rp_upper pr_est pr_lower pr_upper di_est di_lower di_upper"
attr_params=`echo ${rv} ${rp} ${pr} ${DeltaI}`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## WRITE HEADERS & RESULTS TO FILE

header_line=`echo run gmst_file ${eval_headers} ${attr_headers}`
echo $header_line >> $tmp_header_file

model_line=`echo ${model_fnm} ${gmst_fnm} ${eval_params} ${attr_params}`
echo $model_line >> $tmp_result_file

