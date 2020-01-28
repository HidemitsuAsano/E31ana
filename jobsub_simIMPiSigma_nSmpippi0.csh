#!/bin/tcsh -f
set DSTVersion="2"
set Version="4"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simSmpippi0${DSTVersion}/"
set CDSDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"
set CDSDIRSUB="${CDSDIR}simDST_simSmpippi0${DSTVersion}"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simIMPiSigma_nSmpippi0/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simIMPiSigma_nSmpippi0_${Version}"
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

cp hadd_simhist_nSmpippi0.csh $OUTDIRSUB/
cp hadd_sim_nSmpippi0_pippimn.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/simIMpisigma_nSmpippi0_all.root simpost/simIMpisigma_nSmpippi0_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_nSmpippi0_pippimn_all.root simpost/simIMpisigma_nSmpippi0_pippimn_v${Version}.root

@ i = 0
while ($i < 400)   

  set EXEC___="./bin/simIMPiSigma"
  set CONF___="conf/Run78/analyzer_kwsk_sim.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILESP=${DATADIR}"sim_nSmpippi0_0${jobnum}.root"
  
  set CDSFILESP=${CDSDIRSUB}"/simDST_nSmpippi0_0${jobnum}.root"
  
  set OUTFILESP=${OUTDIRSUB}"/simIMpisigma_nSmpippi0_0${jobnum}.root"

  echo ${INPFILESP}
  echo ${CDSFILESP} 
  echo ${OUTFILESP}
  
  set lognamesp = "${logdir}/runnSmpippi0$i.log"
  bsub -o $lognamesp -q s ${EXEC___} ${CONF___} ${OUTFILESP} ${INPFILESP} ${CDSFILESP}
    @ i ++
end

while (1)
  @ njob=`bjobs | wc -l`
  if ( $njob < 1 ) then 
    echo "all jobs finished"
    cd $OUTDIRSUB
    tcsh hadd_simhist_nSmpippi0.csh
    tcsh hadd_sim_nSmpippi0_pippimn.csh
    cd -
    break
  endif
  echo "$njob jobs running" 
  sleep 60
end

echo "aggrigation is finished"
echo "start time ${starttime}"
echo "end time `date '+%y/%m/%d %H:%M:%S'`"
