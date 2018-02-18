#! /bin/sh

vr=1
trackvr=4

#for i in {26..43}
#for i in {44..55}
#for i in {102..240}
#for i in {79..310}
#for i in {79..114}
#for i in {153..163}
#for i in {176..240}
#for i in {253..310} #D2
for i in {1..310}
do
conf=conf/Run62/analyzer_run62_20150610.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/myanabl/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/myanabl "$conf" "$evnum" "$out" "$in"
done
