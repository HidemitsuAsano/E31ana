#!/bin/tcsh -f


root -b -q 'plot_IMLambdaPim.C+'
root -l -b -q 'plot_IMLambdaPim.C+ ("evanaIMLambdaPim_ppimpim_v16.root")'     
root -l -b -q 'plot_IMLambdaPim.C+ ("../simpost/simIMLpim_ppimpim_v20.root")'     
