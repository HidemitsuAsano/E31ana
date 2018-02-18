#! /bin/sh

vr=51
#for H2
trackvr=5
anavr=5
#for D2
#trackvr=4
#anavr=5

#for i in {26..43}
#for i in {44..55}
#for i in {1..310}
#for i in {79..310}
#for i in {79..114}
#for i in {153..163}
#for i in {176..208}
for i in {90..90} #H2
#for i in {90..240} #H2
#for i in {253..310} #D2
#for i in {90..140} #H2
#for i in {141..240} #H2
do
#"conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/MyStudy/analyzer_run62.conf
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
