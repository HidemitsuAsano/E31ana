#!/bin/tcsh -f

@ i = 0
set PARAMDIR="param/Run78/"

while ($i < 812)   

set jobnum=`printf  "%04d"  $i`
#set GainMap="GainMapBL_List:        ${PARAMDIR}GainMapBL/runbyrun_v1/GainMapBL_${jobnum}_xt.param     $i  $i"
#set XTMap="XTMapBL_List:   ${PARAMDIR}XTMapBL/xt_run$i.param     $i  $i"
echo "GainMapBL_List:       ${PARAMDIR}GainMapBL/runbyrun_v1/GainMapBL_${jobnum}_xt.param     $i  $i" 
#echo ${XTMap}
  @ i ++
end 

