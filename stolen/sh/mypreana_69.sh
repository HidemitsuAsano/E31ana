#! /bin/sh

vr=1
trackvr=2
conf=conf/Run69/analyzer.conf

#for i in {1..15}
for i in {120..151}
#for i in {127..127}
#for i in {103..103}
#for i in {149..150}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run69/mypreana/run69_%04d_v%02d.root $i $vr)
in=$(printf e17data/Run69/e62_069_%04d.dat $i)
track=$(printf root/ana/Run69/cdctracking/run69_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreanastop "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done
