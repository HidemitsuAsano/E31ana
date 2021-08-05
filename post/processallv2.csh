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

root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v142.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v142.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v142.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v142.root",3)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v142.root")'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v142.root",1)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v142.root",2)'
root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v142.root",3)'

#root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root")'
#root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root",1)'
#root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root",2)'
#root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSppim_pippimn_v132_MIX_cut4.root",3)'
#root -l -b -q 'plot_IMpisigma.C+ ("../simpost/simIMpisigma_nSmpip_pippimn_v132_MIX_cut4.root")'
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
