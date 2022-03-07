#!/bin/tcsh -f
set Version="1"
set DSTVersion="1"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simpLpimpi0${DSTVersion}/"
set CDSDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"
set CDSDIRSUB="${CDSDIR}simDSTpLpimpi0${DSTVersion}"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simIMLPim_pLpimpi0/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simIMLpim_pLpimpi0_${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}_v${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
else 
 echo "version exist v"${Version}
 exit 0
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/UserSimIMLPim.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/
set CONF___="conf/Run78/analyzer_kwsk_sim.conf"
cp $CONF___ $OUTDIRSUB

cp hadd_simlpimhist.csh $OUTDIRSUB/
cp hadd_simlpim_ppimpim.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/simIMLpim_all.root simpost/simIMLpim_pLpimpi0_v${Version}.root
ln -s $OUTDIRSUB/simIMLpim_ppimpim_all.root simpost/simIMLpim_ppimpim_pLpimpi0_v${Version}.root

@ i = 0
while ($i < 400)   

  set EXEC___="./bin/simIMLPim"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"sim_pLpimpi0_0${jobnum}.root"
  
  set CDSFILE=${CDSDIRSUB}"/simDST_pLpimpi0_0${jobnum}.root"
  
  set OUTFILE=${OUTDIRSUB}"/simIMLpim_0${jobnum}.root"

  echo ${INPFILE}
  echo ${CDSFILE} 
  echo ${OUTFILE}
  
  set logname = "${logdir}/runIMLpim_pLpimpi0$i.log"
  bsub -o $logname -q s ${EXEC___} ${CONF___} ${OUTFILE} ${INPFILE} ${CDSFILE}
  @ i ++
end

while (1)
  @ njob=`bjobs | wc -l`
  if ( $njob < 1 ) then 
    echo "all jobs finished"
    cd $OUTDIRSUB
    tcsh hadd_simlpimhist.csh
    tcsh hadd_simlpim_ppimpim.csh
    cd -
    break
  endif
  echo "$njob jobs running" 
  sleep 60
end

echo "aggrigation is finished"
echo "start time ${starttime}"
echo "end time `date '+%y/%m/%d %H:%M:%S'`"
