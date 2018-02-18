#!/bin/sh
#	root -l -b -q 'merge_check.C ("run47_mhtdc_0093.root", "merge_run47_0093.root")'
i=$1
CONF=conf/Run47/analyzer.conf
FOLDER1=/s/e15/data/Run47
FOLDER2=/w/e15/data/Run47/mtdc_tree
FOLDER3=/w/e15/data/Run47/mtdc
OUTFOLDER=/w/e15/data/Run47/merge_check
PREFIX1=run47_000
PREFIX2=run47_mhtdc_000
PREOUT=merge_run47_000
PROG1=merge_check
PROG2=mtdc_decoder
while [ $i -le $2 ]; do
    export RUN=$i
    if [ $i -ge 10 ]; then
	PREFIX1=run47_00
	PREFIX2=run47_mhtdc_00
	PREOUT=merge_run47_00
	if [ $i -ge 100 ]; then
	    PREFIX1=run47_0
	    PREFIX2=run47_mhtdc_0
	    PREOUT=merge_run47_0
	fi
    fi
    if test -f $FOLDER1/$PREFIX1$i.dat; then
	if test -f $FOLDER2/$PREFIX2$i.root; then
	TKO=$FOLDER1/$PREFIX1$i.dat
	MTDC=$FOLDER2/$PREFIX2$i.root
	nice $PROG1 $CONF $OUTFOLDER/$PREOUT$i.root $TKO $MTDC
	root -l -b -q 'merge_check.C ("'$MTDC'","'$OUTFOLDER/$PREOUT$i.root'")'
	else
	    echo $FOLDER2/$PREFIX2$i.root " not found"
	    echo "start to create mtdc tree ..."
	    nice $PROG2 $FOLDER3/$PREFIX2$i.dat $FOLDER2/$PREFIX2$i.root
	    if test -f $FOLDER2/$PREFIX2$i.root; then
		echo "mtdc tree created."
		echo "start to merge mtdc and tko ..."
		TKO=$FOLDER1/$PREFIX1$i.dat
		MTDC=$FOLDER2/$PREFIX2$i.root
		nice $PROG1 $CONF $OUTFOLDER/$PREOUT$i.root $TKO $MTDC
		root -l -b -q 'merge_check.C ("'$MTDC'","'$OUTFOLDER/$PREOUT$i.root'")'
	    else
		echo "I cannot create mtdc tree file !!"
#		exit
	    fi
	fi
    else
	echo $FOLDER1/$PREFIX1$i.dat " not found"
    fi

    echo $DATA
    i=`expr $i + 1`
done

