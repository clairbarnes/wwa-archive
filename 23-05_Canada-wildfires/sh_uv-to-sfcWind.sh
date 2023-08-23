#!/bin/bash

# script to convert u and v wind vectors to windspeed
# Chaining fails for some reason, so using an intermediate temporary file

fl=`ls crcm5/*/*UU*.nc4`
for u_file in $fl; do
    v_file=${u_file/UU/VV};
    outfile=${u_file/UU/WS};
    file_out=crcm5_wind/${outfile##*/};
    echo $file_out;
    cdo -s chname,uas,sfcWind -sqrt -add -sqr -selname,uas $u_file -sqr -selname,vas $v_file $file_out
done
echo "sfcWind complete"


