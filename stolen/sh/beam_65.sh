#! /bin/sh

vr=3
trackvr=3

#for i in {101..101} #All Data
for i in {1..846} #All Data
#for i in {561..660} #All Data
#for i in {1..208} #All Data
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/beam/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
done
