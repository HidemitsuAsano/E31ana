#!/bin/tcsh -f
set Version="4"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/sim_nS0pipi${Version}/"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simcds"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}simDST_nS0pipi${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/UserSimDatG4.cpp $OUTDIRSUB/

#ln -s $OUTDIRSUB/simIMpisigma_all.root simpost/simIMpisigma_all_v${Version}.root

@ i = 0
while ($i < 400)   

  set EXEC___="./bin/sim"
  set CONF___="conf/Run78/analyzer_kwsk_sim.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"sim_nS0pippim_0${jobnum}.root"
  
  set OUTFILE=${OUTDIRSUB}"/simDST_nS0pippim_0${jobnum}.root"

  echo ${INPFILE}
  echo ${OUTFILE}

  set logname = "${logdir}/run_nS0pippim$i.log"
  bsub -o $logname -q l ${EXEC___} ${CONF___} ${OUTFILE} ${INPFILE} 
    @ i ++
end

