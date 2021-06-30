#!/bin/tcsh -f

hadd -f simIMpisigma_K0nn_all1.root simIMpisigma_K0nn_0[0-9][0-9][0-9].root
hadd -f simIMpisigma_K0nn_all2.root simIMpisigma_K0nn_01[0-9][0-9][0-9].root
hadd -f simIMpisigma_K0nn_all3.root simIMpisigma_K0nn_02[0-9][0-9][0-9].root
hadd -f simIMpisigma_K0nn_all4.root simIMpisigma_K0nn_03[0-9][0-9][0-9].root

hadd -f simIMpisigma_K0nn_all.root simIMpisigma_K0nn_all[1-4].root
