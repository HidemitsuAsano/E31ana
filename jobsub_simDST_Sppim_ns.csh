#!/bin/tcsh -f
set Version="2"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simSppim_ns${Version}/"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simcds"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 
else 
 echo "version exist v"${Version}
 exit 0
endif

set OUTDIRSUB="${OUTDIR}simDSTSppim_ns${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/UserSimDatG4.cpp $OUTDIRSUB/


@ i = 0
while ($i < 400)   

  set EXEC___="./bin/sim"
  set CONF___="conf/Run78/analyzer_kwsk_sim.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"sim_Sppim_ns_0${jobnum}.root"
  
  set OUTFILE=${OUTDIRSUB}"/simDST_Sppim_ns_0${jobnum}.root"

  echo ${INPFILE}
  echo ${OUTFILE}

  set logname = "${logdir}/runSppim_ns$i.log"
  bsub -o $logname -q s ${EXEC___} ${CONF___} ${OUTFILE} ${INPFILE} 
    @ i ++
end

