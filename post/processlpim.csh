#!/bin/tcsh -f

root -b -q 'plot_IMLambdaPim.C+'
root -l -b -q 'plot_IMLambdaPim.C+ ("evanaIMLambdaPim_ppimpim_v25.root")'     
root -l -b -q 'plot_IMLambdaPim.C+ ("../simpost/simIMLpim_ppimpim_v29.root")'     
#root -l -b -q '../simpost/GetAccMapLpim.C'
#root -l -b -q 'CS_IMLambdaPim.C'
#root -l -b -q 'plot_IMLambdaPim.C+ ("../simpost/simIMLpim_ppimpim_pS0pim_v2.root")'     
