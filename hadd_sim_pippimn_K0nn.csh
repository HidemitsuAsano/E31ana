#!/bin/tcsh -f

hadd -f simIMpisigma_K0nn_pippimn_all1.root simIMpisigma_K0nn_0[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_K0nn_pippimn_all2.root simIMpisigma_K0nn_01[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_K0nn_pippimn_all3.root simIMpisigma_K0nn_02[0-9][0-9][0-9]_pippimn.root
hadd -f simIMpisigma_K0nn_pippimn_all4.root simIMpisigma_K0nn_03[0-9][0-9][0-9]_pippimn.root

hadd -f simIMpisigma_K0nn_pippimn_all.root simIMpisigma_K0nn_pippimn_all[1-4].root
