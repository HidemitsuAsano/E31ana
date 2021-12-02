#!/bin/tcsh -f
set Version="18"
set DATADIR="/group/had/knucl/e15/data/Run78/"
set OUTDIR="/group/had/knucl/e15/asano/Run78/"
set KWSKDIR="/group/had/knucl/e15/shinngo/Run78/evtracking/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_IMLambdaPim_${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}IMLambdaPimv${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
else
 echo "version exist v"${Version}
 exit 0 
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/EventAnalysisIMLambdaPim.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/

cp hadd_IMLPimhist.csh $OUTDIRSUB/
cp hadd_IMppimpim.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/evanaIMLambdaPim_all.root post/evanaIMLambdaPim_v${Version}.root
ln -s $OUTDIRSUB/evanaIMLambdaPim_all_ppimpim.root post/evanaIMLambdaPim_ppimpim_v${Version}.root

@ i = 100
while ($i < 813)   

  set EXEC___="./bin/evlpim"
  set CONF___="conf/Run78/analyzer_kwsk.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"run78_0${jobnum}.dat"
  set OUTFILE=${OUTDIRSUB}"/evanaIMLambdaPim_0${jobnum}.root"
  set CDSFILE=${KWSKDIR}"run78_0${jobnum}_evtracking.root"

  echo ${INPFILE}
  echo ${OUTFILE}
  echo ${CDSFILE}

  set logname = "${logdir}/run$i.log"
  bsub -o $logname -q s ${EXEC___} ${CONF___} ${i} ${OUTFILE} ${INPFILE} ${CDSFILE} 
    @ i ++
end

while (1)
  @ njob=`bjobs | wc -l`
  if ( $njob < 1 ) then 
    echo "all jobs finished"
    cd $OUTDIRSUB
    tcsh hadd_IMLPimhist.csh
    tcsh hadd_IMppimpim.csh
    cd - 
    break
  endif
  echo "$njob jobs running" 
  sleep 60
end

echo "aggrigation is finished"
echo "start time ${starttime}"
echo "end time `date '+%y/%m/%d %H:%M:%S'`"
