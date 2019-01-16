#!/bin/tcsh -f
set Version="29"
set DATADIR="/group/had/knucl/e15/data/Run78/"
set OUTDIR="/group/had/knucl/e15/asano/Run78/"
set KWSKDIR="/group/had/knucl/e15/shinngo/Run78/evtracking/"

set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_IMpisigma"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}IMpisigmav${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/EventAnalysisIMPiSigma.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/

cp hadd_IMhist.csh $OUTDIRSUB/
cp hadd_IMnpippim.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/evanaIMpisigma_all.root post/evanaIMpisigma_all_v${Version}.root
ln -s $OUTDIRSUB/evanaIMpisigma_all_npippim.root post/evanaIMpisigma_all_npippim_v${Version}.root

@ i = 100
while ($i < 812)   

  set EXEC___="./bin/evpisigma"
  set CONF___="conf/Run78/analyzer_kwsk.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"run78_0${jobnum}.dat"
  set OUTFILE=${OUTDIRSUB}"/evanaIMpisigma_0${jobnum}.root"
  set CDSFILE=${KWSKDIR}"run78_0${jobnum}_evtracking.root"

  echo ${INPFILE}
  echo ${OUTFILE}
  echo ${CDSFILE}

  set logname = "${logdir}/run$i.log"
  bsub -o $logname -q s ${EXEC___} ${CONF___} ${i} ${OUTFILE} ${INPFILE} ${CDSFILE} 
    @ i ++
end

