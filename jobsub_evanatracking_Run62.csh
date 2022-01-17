#!/bin/tcsh -f

set DATADIR="/group/had/knucl/e15/data/Run62/"
set OUTDIR="/group/had/knucl/e15/asano/Run62/"

if( ! -d $OUTDIR) then 
  mkdir -p  $OUTDIR
endif 

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_cdstrack"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 


@ i = 100
while ($i < 241)   


set EXEC___="./bin/evtracking"
set CONF___="conf/Run78/analyzer_run62_20211019.conf"
set jobnum=`printf  "%03d"  $i`

set INPFILE=${DATADIR}"run62_0${jobnum}.dat"
set OUTFILE=${OUTDIR}"CDStrack_0${jobnum}.root"

echo ${INPFILE}
echo ${OUTFILE}


set logname = "${logdir}/run$i.log"
bsub -o $logname -q l ${EXEC___} ${CONF___} ${i} ${OUTFILE} ${INPFILE} 
  @ i ++
end
