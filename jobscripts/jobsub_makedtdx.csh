#!/bin/tcsh -f

set INDIR="/group/had/knucl/e15/asano/Run78/"

@ i = 0
while ($i < 812)   


set jobnum=`printf  "%03d"  $i`

set INFILE=${INDIR}"evanaxt_0${jobnum}.root"

echo ${INFILE}

root -b -q  'macro/MakeDtDx2_Run78.C('$i')'
  @ i ++
end
