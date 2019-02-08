#!/bin/tcsh -f
set Version="5"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/sim3/"
set CDSDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"
set CDSDIRSUB="${CDSDIR}simDST3"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simIMPiSigma/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simIMPiSigma_DoraAir_${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}DoraAir_v${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/UserSimIMPiSigma.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/

cp hadd_simhist.csh $OUTDIRSUB/
cp hadd_sim_pippimn.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/simIMpisigma_nSppim_all.root simpost/simIMpisigma_nSppim_all_DoraAir_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_nSmpip_all.root simpost/simIMpisigma_nSmpip_all_DoraAir_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_nSppim_pippimn_all.root simpost/simIMpisigma_nSppim_pippimn_all_DoraAir_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_nSmpip_pippimn_all.root simpost/simIMpisigma_nSmpip_pippimn_all_DoraAir_v${Version}.root

@ i = 0
while ($i < 400)   

  set EXEC___="./bin/simIMPiSigma"
  set CONF___="conf/Run78/analyzer_kwsk_sim_DoraAir.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILESP=${DATADIR}"sim_nSppim_0${jobnum}.root"
  set INPFILESM=${DATADIR}"sim_nSmpip_0${jobnum}.root"
  
  set CDSFILESP=${CDSDIRSUB}"/simDST_nSppim_0${jobnum}.root"
  set CDSFILESM=${CDSDIRSUB}"/simDST_nSmpip_0${jobnum}.root"
  
  set OUTFILESP=${OUTDIRSUB}"/simIMpisigma_nSppim_0${jobnum}.root"
  set OUTFILESM=${OUTDIRSUB}"/simIMpisigma_nSmpip_0${jobnum}.root"

  echo ${INPFILESP}
  echo ${CDSFILESP} 
  echo ${OUTFILESP}
  
  set lognamesp = "${logdir}/runSp$i.log"
  bsub -o $lognamesp -q s ${EXEC___} ${CONF___} ${OUTFILESP} ${INPFILESP} ${CDSFILESP}

  echo ${INPFILESM}
  echo ${CDSFILESM} 
  echo ${OUTFILESM}

  set lognamesm = "${logdir}/runSm$i.log"
  bsub -o $lognamesm -q s ${EXEC___} ${CONF___} ${OUTFILESM} ${INPFILESM} ${CDSFILESM}
    @ i ++
end

while (1)
  if ( (`bjobs | wc -l`) < 1 ) then 
    echo "all bjob finished"
    cd $OUTDIRSUB
    tcsh hadd_simhist.csh
    tcsh hadd_sim_pippimn.csh
    break
  endif
  echo "running" 
  sleep 30
end

echo "aggrigation is finished"
