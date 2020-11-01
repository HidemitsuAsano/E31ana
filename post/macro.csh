#!/bin/tcsh -f
set histname="evanaIMpisigma_v193.root"
#echo "${histname}"
root -l -q -b  'plothists.C+("'"${histname}"'")'
