#! /bin/sh

vr=1
trackvr=1

#for i in {1..10}   #BLC2 H.V./Vth scan
#for i in {67..70} #IH Bias scan
#for i in {71..74} #IH Bias scan
#for i in {77..78}  #IH Bias scan
#for i in {80..87}  #CDS Cosmic w/o Dora
#for i in {88..89}  #CDS Cosmic w/ Dora
#for i in {99..99}  #CDS Cosmic w/ Dora
#for i in {2002..2002}  #NC check
#Run62
for i in {1..310}  #H2 run
#for i in {253..310}  #H2 run
do
conf=conf/Run62/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/mypreana/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
done
