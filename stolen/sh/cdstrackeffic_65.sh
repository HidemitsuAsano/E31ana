#! /bin/sh

vr=1
travr=1
conf=conf/Run65/analyzer_run65.conf
for i in {471..471}  #CDC H.V./Vth study
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/cds/run65_%04d_v%02d_trackeffic.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $travr)
bsub -ql -o tmp.log ./bin/cdstrackeffic "$conf" "$evnum" "$out" "$in" "$track"
done
