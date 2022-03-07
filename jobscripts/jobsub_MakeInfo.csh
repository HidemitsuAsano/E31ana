#!/bin/tcsh -f

set DATADIR="/group/had/knucl/e15/data/Run78/"
set OUTDIR="/group/had/knucl/e15/asano/Run78/"
set KWSKDIR="/group/had/knucl/e15/shinngo/Run78/evtracking/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_MakeInfo"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 
@ i = 100
while ($i < 812)   


set EXEC___="./bin/evanaMakeInfo"
set CONF___="conf/Run78/analyzer_kwsk.conf"
set jobnum=`printf  "%03d"  $i`

set INPFILE=${DATADIR}"run78_0${jobnum}.dat"
set OUTFILE=${OUTDIR}"evanaDST_0${jobnum}.root"
set CDSFILE=${KWSKDIR}"run78_0${jobnum}_evtracking.root"

echo ${INPFILE}
echo ${OUTFILE}
echo ${CDSFILE}

set logname = "${logdir}/run$i.log"
bsub -o $logname -q l ${EXEC___} ${CONF___} ${i} ${OUTFILE} ${INPFILE} ${CDSFILE}
  @ i ++
end
