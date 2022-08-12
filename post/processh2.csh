#!/bin/tcsh -f

root -b -q 'plot_IMsigma_h2.C+'
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14.root",2)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",2,0)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",2,1)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",2,-1)'     
root -l -b -q 'SubtractMix_H2.C (2,0)'
root -l -b -q 'SubtractMix_H2.C (2,1)'
root -l -b -q 'SubtractMix_H2.C (2,-1)'
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14.root",4)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",4,0)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",4,1)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",4,-1)'     
root -l -b -q 'SubtractMix_H2.C (4,0)'
root -l -b -q 'SubtractMix_H2.C (4,1)'
root -l -b -q 'SubtractMix_H2.C (4,-1)'
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14.root",6)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",6,0)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",6,1)'     
root -l -b -q 'plot_IMsigma_h2.C+ ("evanaIMsigma_npi_h2_v14_MIX.root",6,-1)'     
root -l -b -q 'SubtractMix_H2.C (6)'
root -l -b -q 'SubtractMix_H2.C (6,1)'
root -l -b -q 'SubtractMix_H2.C (6,-1)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_Sppim_npi_v15.root",2)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_Smpip_npi_v15.root",2)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_K0n_npi_v2.root",2)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_pipiL_npi_v1.root",2)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_Sppim_npi_v15.root",4)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_Smpip_npi_v15.root",4)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_K0n_npi_v2.root",4)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_pipiL_npi_v1.root",4)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_Sppim_npi_v15.root",6)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_Smpip_npi_v15.root",6)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_K0n_npi_v2.root",6)'
root -l -b -q 'plot_IMsigma_h2.C+ ("../simpost/simIMsigma_H2_pipiL_npi_v1.root",6)'
root -l -b -q 'GetAccH2.C (15,2)'
root -l -b -q 'GetAccH2.C (15,4)'
root -l -b -q 'GetAccH2.C (15,6)'
root -l -b -q 'CS_sigma_h2.C'
root -l       'comp_pastDataH2.C'