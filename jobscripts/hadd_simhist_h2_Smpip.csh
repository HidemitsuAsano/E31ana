#!/bin/tcsh -f

hadd -f simIMsigma_H2_Smpip_all1.root simIMsigma_H2_Smpip_0[0-9][0-9][0-9].root
hadd -f simIMsigma_H2_Smpip_all2.root simIMsigma_H2_Smpip_01[0-9][0-9][0-9].root

hadd -f simIMsigma_H2_Smpip_all.root simIMsigma_H2_Smpip_all[1-2].root 
