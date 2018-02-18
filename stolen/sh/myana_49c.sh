#! /bin/sh

#Run49c
vr=2
anavr=10
for i in {60..156} #H2 production
do
conf=conf/Run49c/analyzer_run49c.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run49c/myana/run49c_%04d_v%02d.root $i $vr)
in=$(printf data/Run49c/run49c_%04d.dat $i)
all=$(printf root/ana/Run49c/myall/run49c_%04d_v%02d.root $i $anavr)
track=$(printf /group/had/knucl/e15/hashimoto/run49c/tra_20131120_%d.root $i)
bsub -qe -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
