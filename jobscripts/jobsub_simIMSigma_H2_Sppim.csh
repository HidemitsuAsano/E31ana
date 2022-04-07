#!/bin/tcsh -f
set Version="15"
set DSTVersion="2"
set DATADIR="/gpfs/group/had/knucl/e15/asano/sim/simSppim${DSTVersion}/"
set CDSDIR="/gpfs/group/had/knucl/e15/asano/sim/simcds/"
set CDSDIRSUB="${CDSDIR}simDSTSppim${DSTVersion}"
set OUTDIR="/gpfs/group/had/knucl/e15/asano/sim/simIMSigma_H2/"

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

cd ../
set SRCDIR="/gpfs/home/had/hiasano/ana/k18ana/src/"
set CSHDIR="/gpfs/home/had/hiasano/ana/k18ana/jobscripts/"
cp $SRCDIR/UserSimIMSigma_H2.cpp $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaAnaPar.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaHist.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.h $OUTDIRSUB/
cp $SRCDIR/IMPiSigmaUtil.cpp $OUTDIRSUB/

cp $CSHDIR/hadd_simhist_h2_Sppim.csh $OUTDIRSUB/
cp $CSHDIR/hadd_sim_h2_Sppim_npi.csh $OUTDIRSUB/
ln -s $OUTDIRSUB/simIMsigma_H2_Sppim_all.root simpost/simIMsigma_H2_Sppim_v${Version}.root
ln -s $OUTDIRSUB/simIMsigma_H2_Sppim_npi_all.root simpost/simIMsigma_H2_Sppim_npi_v${Version}.root

@ i = 0
while ($i < 1600)   

  set EXEC___="./bin/simIMSigma_h2"
  set CONF___="conf/Run78/analyzer_run62_20211019.conf"
  set jobnum=`printf  "%03d"  $i`

  set INPFILESP=${DATADIR}"sim_Sppim_0${jobnum}.root"
  
  set CDSFILESP=${CDSDIRSUB}"/simDST_Sppim_0${jobnum}.root"
  
  set OUTFILESP=${OUTDIRSUB}"/simIMsigma_H2_Sppim_0${jobnum}.root"

  echo ${INPFILESP}
  echo ${CDSFILESP} 
  echo ${OUTFILESP}
  
  set lognamesp = "${logdir}/runSp$i.log"
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
    tcsh hadd_sim_h2_Sppim_npi.csh
    tcsh hadd_simhist_h2_Sppim.csh
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
