#! /bin/sh

vr=3
trackvr=2
conf=conf/Run69/analyzer.conf

for i in {100..167} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run69/beam/run69_%04d_v%02d.root $i $vr)
in=$(printf e17data/Run69/e62_069_%04d.dat $i)
bsub -qs -o tmp.log ./bin/beamstop "$conf" "$evnum" "$out" "$in"
done
