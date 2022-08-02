#!/bin/tcsh -f

root -b -q 'plot_IMpisigma.C+'
set VERSION=241

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v241.root",3,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v241_MIX_cut4.root",3,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v241_MIX_cut4.root",3,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v241_MIX_cut4.root",3,2,1)' 

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v156.root",3,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v156.root",3,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",3,2)'

set hname=`hostname -s`
echo ${hname}
while (1)
  @ njob=`bjobs | grep ${hname} | wc -l`
  if ( $njob < 1 ) then 
    root -l -b -q 'SubtractMix.C'
    break
  endif
  echo "$njob jobs running" 
  sleep 30
end

cd ../simpost/
root -l -b -q 'GetAccMap.C (2)'
cd ../post/
root -l -b -q 'plot_AfterDecompos.C(2,0)'
root -l -b -q 'plot_AfterDecompos.C(2,1)'
root -l -b -q 'plot_AfterDecompos.C(2,-1)'

root -l -b -q 'CS_finals.C'
#mv *.pdf pdf/
