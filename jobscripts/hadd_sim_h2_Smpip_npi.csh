#!/bin/tcsh -f

hadd -f simIMsigma_H2_Smpip_npi_all1.root simIMsigma_H2_Smpip_0[0-9][0-9][0-9]_npi.root
hadd -f simIMsigma_H2_Smpip_npi_all2.root simIMsigma_H2_Smpip_01[0-9][0-9][0-9]_npi.root

hadd -f simIMsigma_H2_Smpip_npi_all.root simIMsigma_H2_Smpip_npi_all[1-2].root 
