#!/bin/sh
#./bin/evpisigma conf/Run78/analyzer_kwsk.conf 601 test_0601.root /group/had/knucl/e15/data/Run78/run78_0601.dat /group/had/knucl/e15/shinngo/Run78/evtracking/run78_0601_evtracking.root  
./bin/evanaraw conf/Run78/analyzer_kwsk.conf 601 test_0601_raw.root run78_0601.dat run78_0601_evtracking.root  
#./bin/evpisigma conf/Run78/analyzer_kwsk_CDHgainv32.conf 601 test_0601.root /group/had/knucl/e15/data/Run78/run78_0601.dat /group/had/knucl/e15/shinngo/Run78/evtracking/run78_0601_evtracking.root  

