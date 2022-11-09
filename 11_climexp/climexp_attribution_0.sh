
model_list=`echo p_lake-chad_jjas-mean_IPSL-CM5A-LR_r1_REMO2009:tas_Amon_IPSL-CM5A-LR_rcp85_r1i1p1 p_lake-chad_jjas-mean_MPI-ESM-LR_r1_REMO2009:tas_Amon_MPI-ESM-LR_rcp85_r1i1p1 p_lake-chad_jjas-mean_MPI-ESM-MR_r1_RegCM4-3:tas_Amon_MPI-ESM-MR_rcp85_r1i1p1 p_lake-chad_jjas-mean_HadGEM2-ES_r1_RegCM4-3:tas_Amon_HadGEM2-ES_rcp85_r1i1p1 p_lake-chad_jjas-mean_HadGEM2-ES_r1_REMO2015:tas_Amon_HadGEM2-ES_rcp85_r1i1p1 p__lake-chad_jjas-meanHadGEM2-ES_r1_RegCM4-7:tas_Amon_HadGEM2-ES_rcp85_r1i1p1 p_lake-chad_jjas-mean_MIROC5_r1_REMO2009:tas_Amon_MIROC5_rcp85_r1i1p1 p_lake-chad_jjas-mean_HadGEM2-ES_r1_REMO2009:tas_Amon_HadGEM2-ES_rcp85_r1i1p1 p_lake-chad_jjas-mean_EC-EARTH_r12_REMO2009:tas_Amon_EC-EARTH_rcp85_r12i1p1 p_lake-chad_jjas-mean_NorESM1-M_r1_RegCM4-7:tas_Amon_NorESM1-M_rcp85_r1i1p1 p_lake-chad_jjas-mean_MPI-ESM-MR_r1_RegCM4-7:tas_Amon_MPI-ESM-MR_rcp85_r1i1p1`

results_file=results_lake-chad.txt
attr_script=climexp_attribution.sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## DOWNLOAD LIST OF AVAILABLE DATA

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
rm tmp_headers.txt tmp_results.txt model_eval.log model_attr.log model_proj.log climexp_uploads.txt