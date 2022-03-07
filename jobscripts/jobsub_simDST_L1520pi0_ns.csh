#!/bin/tcsh -f
set Version="1"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simL1520pi0_ns${Version}/"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simcds"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}simDSTL1520pi0_ns${Version}"
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

  set INPFILE=${DATADIR}"sim_L1520pi0_ns_0${jobnum}.root"
  
  set OUTFILE=${OUTDIRSUB}"/simDST_L1520pi0_ns_0${jobnum}.root"

  echo ${INPFILE}
  echo ${OUTFILE}

  set logname = "${logdir}/runL1520pi0$i.log"
  bsub -o $logname -q s ${EXEC___} ${CONF___} ${OUTFILE} ${INPFILE} 
    @ i ++
end

