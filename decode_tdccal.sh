#!/bin/sh

#taken in Jan.15 ,2018 before starting Run78 beam time 
for i in {01..40}
do 
echo $i 
$ANAHOME/bin/evanatko conf/Run78/analyzer.conf $i tkoout/out$i.root $HOME/Run78predata/tcal_00$i.dat
done

#taken in Jan.26 , 2018 
$ANAHOME/bin/evanatko conf/Run78/analyzer.conf $i tkoout/out163.root $HOME/Run78data/run78_0163.dat
$ANAHOME/bin/evanatko conf/Run78/analyzer.conf $i tkoout/out164.root $HOME/Run78data/run78_0164.dat
$ANAHOME/bin/evanatko conf/Run78/analyzer.conf $i tkoout/out165.root $HOME/Run78data/run78_0165.dat


#taken in Feb. 13th, 2018
for i in {01..42}
do 
echo $i 
$ANAHOME/bin/evanatko conf/Run78/analyzer.conf $i tkoout/out$i_2.root $HOME/Run78middata/tcal_00$i.dat
done
