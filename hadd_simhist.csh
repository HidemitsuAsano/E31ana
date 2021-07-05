#!/bin/tcsh -f

hadd -f simIMpisigma_nSmpip_all1.root simIMpisigma_nSmpip_0[0-9][0-9][0-9].root
hadd -f simIMpisigma_nSppim_all1.root simIMpisigma_nSppim_0[0-9][0-9][0-9].root
hadd -f simIMpisigma_nSmpip_all2.root simIMpisigma_nSmpip_01[0-9][0-9][0-9].root
hadd -f simIMpisigma_nSppim_all2.root simIMpisigma_nSppim_01[0-9][0-9][0-9].root
hadd -f simIMpisigma_nSmpip_all3.root simIMpisigma_nSmpip_02[0-9][0-9][0-9].root
hadd -f simIMpisigma_nSppim_all3.root simIMpisigma_nSppim_02[0-9][0-9][0-9].root
hadd -f simIMpisigma_nSmpip_all4.root simIMpisigma_nSmpip_03[0-9][0-9][0-9].root
hadd -f simIMpisigma_nSppim_all4.root simIMpisigma_nSppim_03[0-9][0-9][0-9].root

hadd -f simIMpisigma_nSppim_all.root simIMpisigma_nSppim_all[1-4].root 
hadd -f simIMpisigma_nSmpip_all.root simIMpisigma_nSmpip_all[1-4].root 
