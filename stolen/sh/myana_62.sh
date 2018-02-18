#! /bin/sh

conf=conf/Run62/analyzer_run62_H2.conf

#Run62
vr=1
anavr=10
trackvr=5
for i in {90..240} #H2 production
#trackvr=4
#for i in {241..310} #D2 production
do
#conf=conf/Run62/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/myana/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
all=$(printf root/ana/Run62/myall/run62_%04d_v%02d.root $i $anavr)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
bsub -ql -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
