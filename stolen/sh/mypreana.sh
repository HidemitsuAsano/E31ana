#! /bin/sh

vr=11
trackvr=1

#for i in {1..10}   #BLC2 H.V./Vth scan
#for i in {67..70} #IH Bias scan
#for i in {71..74} #IH Bias scan
#for i in {77..78}  #IH Bias scan
#for i in {80..87}  #CDS Cosmic w/o Dora
#for i in {88..89}  #CDS Cosmic w/ Dora
#for i in {99..99}  #CDS Cosmic w/ Dora
#for i in {2002..2002}  #NC check
#Run64
for i in {1..80}  #pi scat.
#for i in {31..33}  #pi scat.
#for i in {50..58}  #pi scat.
#for i in {80..80}  #pi scat.
#Run62
#for i in {153..163}  #H2 run
do
conf=conf/Run64/analyzer_run64.conf
#conf=conf/Run62/analyzer.hashi.20150723
evnum=$(printf %04d $i)
out=$(printf root/ana/Run64/mypreana/run64_%04d_v%02d.root $i $vr)
in=$(printf data/Run64/run64_%04d.dat $i)
track=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
done
