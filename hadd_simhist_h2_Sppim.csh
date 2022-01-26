#!/bin/tcsh -f

hadd -f simIMsigma_H2_Sppim_all1.root simIMsigma_H2_Sppim_0[0-9][0-9][0-9].root
hadd -f simIMsigma_H2_Sppim_all2.root simIMsigma_H2_Sppim_01[0-9][0-9][0-9].root
hadd -f simIMsigma_H2_Sppim_all3.root simIMsigma_H2_Sppim_02[0-9][0-9][0-9].root
hadd -f simIMsigma_H2_Sppim_all4.root simIMsigma_H2_Sppim_03[0-9][0-9][0-9].root

hadd -f simIMsigma_H2_Sppim_all.root simIMsigma_H2_Sppim_all[1-4].root 
