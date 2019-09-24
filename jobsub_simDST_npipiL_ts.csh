#!/bin/tcsh -f
set Version="1"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simnpipiL_ts${Version}/"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simcds"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}simDSTnpipiL_ts${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/UserSimDatG4.cpp $OUTDIRSUB/

#ln -s $OUTDIRSUB/simIMpisigma_all.root simpost/simIMpisigma_all_v${Version}.root

@ i = 0
while ($i < 400)   

  set EXEC___="./bin/sim"
  set CONF___="conf/Run78/analyzer_kwsk_sim_DoraAir.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"sim_npipiL_ts_0${jobnum}.root"
  
  set OUTFILE=${OUTDIRSUB}"/simDST_npipiL_ts_0${jobnum}.root"

  echo ${INPFILE}
  echo ${OUTFILE}

  set lognamesp = "${logdir}/runnpipiL_ts$i.log"
  bsub -o $lognamesp -q l ${EXEC___} ${CONF___} ${OUTFILE} ${INPFILE} 
    @ i ++
end

