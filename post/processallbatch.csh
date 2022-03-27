#!/bin/tcsh -f

root -b -q 'plot_IMpisigma.C+'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",0,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",1,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",2,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",3,2)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,2,0)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,2,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,2,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,2,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,2,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,2,1)' 

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",0,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",1,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",2,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",3,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",0,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",1,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",2,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",3,2)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",0,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",1,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",2,2)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",3,2)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",0,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",1,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",2,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",3,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,4)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,4,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,4,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,4,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,4,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,4,1)' 

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",0,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",1,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",2,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",3,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",0,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",1,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",2,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",3,4)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",0,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",1,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",2,4)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",3,4)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",0,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",1,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",2,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239.root",3,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,6)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,6,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,6,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,6,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,6,-1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",0,6,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",1,6,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",2,6,1)' 
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v239_MIX_cut4.root",3,6,1)' 

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",0,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",1,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",2,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v152.root",3,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",0,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",1,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",2,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v152.root",3,6)'

bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",0,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",1,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",2,6)'
bsub -q s root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v27.root",3,6)'


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
#cd ../simpost/
#mv *.pdf pdf/
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_npipiL_pippimn_v23.root")' 
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_pipiL_ns_pippimn_v1.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_npipiL_ts_pippimn_v1.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nS0pippim_pippimn_v5.root")' 
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_S0pippim_ns_pippimn_v1.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppimpi0_pippimn_v4.root")' 
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpippi0_pippimn_v4.root")' 
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_Sppim_ns_pippimn_v2.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_Smpip_ns_pippimn_v2.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0n_ns_pippimn_v12.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v8.root")' 
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0_nnts_pippimn_v10.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0_Knts_pippimn_v3.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0_Kmpts_pippimn_v1.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nL1520_pippimn_v2.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_L1520pi0_ns_pippimn_v1.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_pipS1385m_ns_pippimn_v1.root")'
#  root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_pimS1385p_ns_pippimn_v1.root")'
