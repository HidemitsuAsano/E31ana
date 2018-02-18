#! /bin/sh

vr=10
trackvr=1

for i in {60..156}
do
conf=conf/Run49c/anahashi20140527.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run49c/mypreana/run49c_%04d_v%02d.root $i $vr)
in=$(printf data/Run49c/run49c_%04d.dat $i)
track=$(printf root/ana/Run49c/cdctracking/run49c_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
done
