#! /bin/sh

#Run49c
vr=1
for i in {60..156} #3He production
do
conf=conf/Run49c/analyzer_run49c.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run49c/myall/run49c_%04d_v%02d.root $i $vr)
in=$(printf data/Run49c/run49c_%04d.dat $i)
track=$(printf /group/had/knucl/e15/hashimoto/run49c/tra_20131120_%d.root $i)
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
