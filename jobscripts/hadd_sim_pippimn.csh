#!/bin/tcsh -f

hadd -f simIMpisigma_nSmpip_pippimn_all1.root simIMpisigma_nSmpip_0[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_nSppim_pippimn_all1.root simIMpisigma_nSppim_0[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_nSmpip_pippimn_all2.root simIMpisigma_nSmpip_01[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_nSppim_pippimn_all2.root simIMpisigma_nSppim_01[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_nSmpip_pippimn_all3.root simIMpisigma_nSmpip_02[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_nSppim_pippimn_all3.root simIMpisigma_nSppim_02[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_nSmpip_pippimn_all4.root simIMpisigma_nSmpip_03[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_nSppim_pippimn_all4.root simIMpisigma_nSppim_03[0-9][0-9][0-9]_pippimn.root

hadd -f simIMpisigma_nSmpip_pippimn_all.root simIMpisigma_nSmpip_pippimn_all[1-4].root 
hadd -f simIMpisigma_nSppim_pippimn_all.root simIMpisigma_nSppim_pippimn_all[1-4].root 
