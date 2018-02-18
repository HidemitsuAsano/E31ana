#! /bin/sh

vr=9
#for H2
trackvr=5
anavr=5

for i in {90..240} #H2
do
#for H2
conf=conf/Run62/analyzer.hashi.20150723
#for D2
#conf=conf/Run62/calib.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/myana/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
ana=$(printf root/ana/Run62/myall/run62_%04d_v%02d.root $i $anavr)
bsub -qe -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$ana"
done
