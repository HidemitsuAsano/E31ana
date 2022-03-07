#!/bin/tcsh -f
set Version="1"
set DSTVersion="1"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simpipiL${DSTVersion}/"
set CDSDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"
set CDSDIRSUB="${CDSDIR}simDSTpipiL${DSTVersion}"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simIMSigma_H2_pipiL/"

set starttime=`date '+%y/%m/%d %H:%M:%S'`
set logbasedir="/home/had/hiasano/logs/"
set date=`date +%Y%m%d_%H%M`
set logdir="${logbasedir}${date}_simIMSigma_H2${Version}"
echo "log files  ${logdir}"

if( ! -d $logdir) then 
  mkdir -p  $logdir
endif 

set OUTDIRSUB="${OUTDIR}v${Version}"
if( ! -d $OUTDIRSUB) then 
 mkdir -p $OUTDIRSUB
endif

set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
cp $SRCDIR/UserSimIMSigma_H2.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/

cp hadd_simhist_h2_pipiL.csh $OUTDIRSUB/
cp hadd_sim_h2_pipiL_npi.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/simIMsigma_H2_pipiL_all.root simpost/simIMsigma_H2_pipiL_v${Version}.root
ln -s $OUTDIRSUB/simIMsigma_H2_pipiL_npi_all.root simpost/simIMsigma_H2_pipiL_npi_v${Version}.root

@ i = 0
while ($i < 400)   

  set EXEC___="./bin/simIMSigma_h2"
  set CONF___="conf/Run78/analyzer_run62_20211019.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILESP=${DATADIR}"sim_pipiL_0${jobnum}.root"
  
  set CDSFILESP=${CDSDIRSUB}"/simDST_pipiL_0${jobnum}.root"
  
  set OUTFILESP=${OUTDIRSUB}"/simIMsigma_H2_pipiL_0${jobnum}.root"

  echo ${INPFILESP}
  echo ${CDSFILESP} 
  echo ${OUTFILESP}
  
  set lognamesp = "${logdir}/runK0$i.log"
  bsub -o $lognamesp -q s ${EXEC___} ${CONF___} ${OUTFILESP} ${INPFILESP} ${CDSFILESP}
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
    tcsh hadd_sim_h2_pipiL_npi.csh
    tcsh hadd_simhist_h2_pipiL.csh
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
