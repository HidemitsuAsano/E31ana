#!/bin/tcsh -f
set Version="23"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/sim${Version}/"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simcds"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}simDST${Version}"
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

  set INPFILESP=${DATADIR}"sim_nSppim_0${jobnum}.root"
  set INPFILESM=${DATADIR}"sim_nSmpip_0${jobnum}.root"
  
  set OUTFILESP=${OUTDIRSUB}"/simDST_nSppim_0${jobnum}.root"
  set OUTFILESM=${OUTDIRSUB}"/simDST_nSmpip_0${jobnum}.root"

  echo ${INPFILESP}
  echo ${OUTFILESP}
  echo ${INPFILESM}
  echo ${OUTFILESM}

  set lognamesp = "${logdir}/runSp$i.log"
  bsub -o $lognamesp -q s ${EXEC___} ${CONF___} ${OUTFILESP} ${INPFILESP} 
  set lognamesm = "${logdir}/runSm$i.log"
  bsub -o $lognamesm -q s ${EXEC___} ${CONF___} ${OUTFILESM} ${INPFILESM} 
    @ i ++
end

