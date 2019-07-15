#!/bin/tcsh -f
set Version="1"
set DSTVersion="1"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simSppim_ns${DSTVersion}/"
set CDSDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"
set CDSDIRSUB="${CDSDIR}simDSTSppim_ns${DSTVersion}"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simIMPiSigma_Sppim_ns/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simIMPiSigma_Sppim_ns_DoraAir_${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}DoraAir_v${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
else 
 echo "version exist v"${Version}
 exit 0
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/UserSimIMPiSigma.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/
set CONF___="conf/Run78/analyzer_kwsk_sim_DoraAir.conf"
cp $CONF___ $OUTDIRSUB

cp hadd_simhist_Sppim_ns.csh $OUTDIRSUB/
cp hadd_sim_pippimn_Sppim_ns.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/simIMpisigma_Sppim_ns_all.root simpost/simIMpisigma_Sppim_ns_DoraAir_v${Version}.root
ln -s $OUTDIRSUB/simIMpisigma_Sppim_ns_pippimn_all.root simpost/simIMpisigma_Sppim_ns_pippimn_DoraAir_v${Version}.root

@ i = 0
while ($i < 400)   

  set EXEC___="./bin/simIMPiSigma"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"sim_Sppim_ns_0${jobnum}.root"
  
  set CDSFILE=${CDSDIRSUB}"/simDST_Sppim_ns_0${jobnum}.root"
  
  set OUTFILE=${OUTDIRSUB}"/simIMpisigma_Sppim_ns_0${jobnum}.root"

  echo ${INPFILE}
  echo ${CDSFILE} 
  echo ${OUTFILE}
  
  set lognamesp = "${logdir}/runSppim_ns$i.log"
  bsub -o $lognamesp -q s ${EXEC___} ${CONF___} ${OUTFILE} ${INPFILE} ${CDSFILE}

  @ i ++
end

while (1)
  @ njob=`bjobs | wc -l`
  if ( $njob < 1 ) then 
    echo "all jobs finished"
    cd $OUTDIRSUB
    tcsh hadd_simhist_Sppim_ns.csh
    tcsh hadd_sim_pippimn_Sppim_ns.csh
    cd -
    break
  endif
  echo "$njob jobs running" 
  sleep 60
end

echo "aggrigation is finished"
echo "start time ${starttime}"
echo "end time `date '+%y/%m/%d %H:%M:%S'`"
