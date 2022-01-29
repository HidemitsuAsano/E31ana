#!/bin/tcsh -f
set Version="146"
set DSTVersion="29"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/sim${DSTVersion}/"
set CDSDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"
set CDSDIRSUB="${CDSDIR}simDST${DSTVersion}"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simIMPiSigma/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simIMPiSigma_${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}v${Version}"
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
ln -s $OUTDIRSUB/simIMpisigma_nSppim_all.root simpost/simIMpisigma_nSppim_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_nSmpip_all.root simpost/simIMpisigma_nSmpip_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_nSppim_pippimn_all.root simpost/simIMpisigma_nSppim_pippimn_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_nSmpip_pippimn_all.root simpost/simIMpisigma_nSmpip_pippimn_v${Version}.root

@ i = 0
while ($i < 3200)   

  set EXEC___="./bin/simIMPiSigma"
  set CONF___="conf/Run78/analyzer_kwsk_sim.conf"
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

set hname=`hostname -s`
echo ${hname}
while (1)
  @ njob=`bjobs | grep ${hname} | wc -l`
  if ( $njob < 1 ) then 
    set agstarttime=`date '+%y/%m/%d %H:%M:%S'`
    echo "all jobs finished"
    cd $OUTDIRSUB
    tcsh hadd_sim_pippimn.csh
    tcsh hadd_simhist.csh
    cd -
    break
  endif
  echo "$njob jobs running" 
  sleep 60
end

echo "aggrigation is finished"
echo "start time ${starttime}"
echo "ag. start time  ${agstarttime}"
echo "end time `date '+%y/%m/%d %H:%M:%S'`"
