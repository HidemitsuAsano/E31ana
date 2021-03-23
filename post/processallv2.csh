#!/bin/tcsh -f

root -b -q 'plot_IMpisigma.C+'
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202.root")' 
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202.root",1)' 
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202.root",2)' 
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202.root",3)' 
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202_MIX_cut4.root")' 
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202_MIX_cut4.root",1)' 
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202_MIX_cut4.root",2)' 
#root -l -b -q 'plot_IMpisigma.C+ ("evanaIMpisigma_npippim_v202_MIX_cut4.root",3)' 

root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132.root",3)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132.root",3)'

root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root",3)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132_MIX_cut4.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132_MIX_cut4.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132_MIX_cut4.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132_MIX_cut4.root",3)'

root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11.root",3)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11_MIX_cut4.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11_MIX_cut4.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11_MIX_cut4.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_K0nn_pippimn_v11_MIX_cut4.root",3)'

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
