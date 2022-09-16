#!/bin/tcsh -f

root -b -q 'plot_IMpisigma.C+'
set VERSION=245

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",0,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",1,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",2,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",3,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",0,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",1,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",2,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",3,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",0,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",1,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",2,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",3,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",0,2,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",1,2,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",2,2,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",3,2,1)' 

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",0,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",1,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",2,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",3,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",0,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",1,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",2,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",3,2)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",0,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",1,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",2,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",3,2)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",0,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",1,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",2,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245.root",3,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",0,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",1,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",2,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",3,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",0,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",1,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",2,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",3,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",0,4,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",1,4,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",2,4,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v245_MIX_cut4.root",3,4,1)' 

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",0,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",1,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",2,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v158.root",3,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",0,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",1,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",2,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v158.root",3,4)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",0,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",1,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",2,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v30.root",3,4)'



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

#mv *.pdf pdf/
