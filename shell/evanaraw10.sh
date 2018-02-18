#!/bin/sh

DIR1="/group/had/knucl/e15/data/Run78/"
DIR2="/group/had/knucl/e15/shinngo/Run78/"
DIR3="/group/had/knucl/e15/shinngo/Run74/evtracking/"

for i in {92..99}
do

for j in {0..99}
do


EXEC___="./bin/evanaraw"
CONF___="conf/Run78/analyzer_"${j}".conf"
INPFILE=${DIR1}"run78_00"${i}".dat"
OUTFILE=${DIR2}"run78_00"${i}"_"${j}".root"
CDCFILE=${DIR3}"run74_0034_evtracking.root"

bsub -q s ${EXEC___} ${CONF___} ${i} ${OUTFILE} ${INPFILE} ${CDCFILE}
#sleep 90s

done
done


