#!/bin/tcsh -f
set Version="245"
set DATADIR="/group/had/knucl/e15/data/Run78/"
set OUTDIR="/group/had/knucl/e15/asano/Run78/"
set KWSKDIR="/group/had/knucl/e15/shinngo/Run78/evtracking/v9/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_IMpisigma_${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}IMpisigmav${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
else 
 echo "version exist v"${Version}
 exit 0
endif

cd ../
set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
set CSHDIR="/gpfs/home/had/hiasano/ana/k18ana/jobscripts/"
cp $SRCDIR/EventAnalysisIMPiSigma.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/
cp conf/Run78/analyzer_kwsk.conf $OUTDIRSUB/

cp $CSHDIR/hadd_IMhist.csh $OUTDIRSUB/
cp $CSHDIR/hadd_IMnpippim.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/evanaIMpisigma_all.root post/evanaIMpisigma_v${Version}.root
ln -s $OUTDIRSUB/evanaIMpisigma_all_npippim.root post/evanaIMpisigma_npippim_v${Version}.root

@ i = 100
while ($i < 813)   

  set EXEC___="./bin/evpisigma"
  set CONF___="./conf/Run78/analyzer_kwsk.conf"
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

set hname=`hostname -s`
echo ${hname}
while (1)
  @ njob=`bjobs | grep ${hname} | wc -l`
  if ( $njob < 1 ) then 
    echo "all jobs finished"
    cd $OUTDIRSUB
    tcsh hadd_IMnpippim.csh
    tcsh hadd_IMhist.csh
    cd - 
    break
  endif
  echo "$njob jobs running" 
  sleep 60
end

cd post/
set histname = "evanaIMpisigma_v${Version}.root"
root -l -q -b  'plothists.C("'"${histname}"'")'
echo "aggrigation is finished"
echo "start time ${starttime}"
echo "end time `date '+%y/%m/%d %H:%M:%S'`"
