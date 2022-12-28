#!/bin/tcsh -f

cd ../simpost/
pwd
root -l -b -q 'GetAccMap.C(2)'
#root -l -b -q 'GetAccMap.C(4)'
#root -l -b -q 'GetAccMap.C(6)'
cd ../post/
root -l -b -q 'plot_AfterDecompos.C(2,0)'
root -l -b -q 'plot_AfterDecompos.C(2,1)'
root -l -b -q 'plot_AfterDecompos.C(2,-1)'
#root -l -b -q 'plot_AfterDecompos.C(4,0)'
#root -l -b -q 'plot_AfterDecompos.C(6,0)'


root -l -b -q 'CS_finals.C'
