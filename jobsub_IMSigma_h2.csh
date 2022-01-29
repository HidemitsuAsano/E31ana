#!/bin/tcsh -f
set Version="4"
set DATADIR="/group/had/knucl/e15/data/Run62/"
set OUTDIR="/group/had/knucl/e15/asano/Run62/"
set CDSTRACKDIR="/group/had/knucl/e15/asano/Run62/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_sigma_h2_${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}IMsigma_h2_v${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
else 
 echo "version exist v"${Version}
 exit 0
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/EventAnalysisIMSigma_H2.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/
cp conf/Run78/analyzer_run62_20211019.conf $OUTDIRSUB/

cp hadd_IMhist_h2.csh $OUTDIRSUB/
cp hadd_IMnpi_h2.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/evanaIMsigma_all.root post/evanaIMsigma_h2_v${Version}.root
ln -s $OUTDIRSUB/evanaIMsigma_all_npi.root post/evanaIMsigma_npi_h2_v${Version}.root

@ i = 100
while ($i < 241)   

  set EXEC___="./bin/evsigma_h2"
  set CONF___="conf/Run78/analyzer_run62_20211019.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILE=${DATADIR}"run62_0${jobnum}.dat"
  set OUTFILE=${OUTDIRSUB}"/evanaIMsigma_h2_0${jobnum}.root"
  set CDSFILE=${CDSTRACKDIR}"CDStrack_0${jobnum}.root"

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
    tcsh hadd_IMnpi_h2.csh
    tcsh hadd_IMhist_h2.csh
    cd - 
    break
  endif
  echo "$njob jobs running" 
  sleep 60
end

#cd post/
#set histname = "evanaIMpisigma_v${Version}.root"
#root -l -q -b  'plothists.C("'"${histname}"'")'
#cd - 
echo "aggrigation is finished"
echo "start time ${starttime}"
echo "end time `date '+%y/%m/%d %H:%M:%S'`"
