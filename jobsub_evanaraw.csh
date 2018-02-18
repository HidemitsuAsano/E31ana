#!/bin/tcsh -f

set DATADIR="/group/had/knucl/e15/data/Run78/"
set OUTDIR="/group/had/knucl/e15/asano/Run78/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_raw"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 
@ i = 0
while ($i < 600)   


set EXEC___="./bin/evanaraw"
set CONF___="conf/Run78/analyzer.conf"
set jobnum=`printf  "%03d"  $i`

set INPFILE=${DATADIR}"run78_0${jobnum}.dat"
set OUTFILE=${OUTDIR}"evanaraw_0${jobnum}.root"

echo ${INPFILE}
echo ${OUTFILE}


set logname = "${logdir}/run$i.log"
bsub -o $logname -q s ${EXEC___} ${CONF___} ${i} ${OUTFILE} ${INPFILE} 
  @ i ++
end
