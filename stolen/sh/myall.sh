#! /bin/sh

#for H2
#vr=5
#trackvr=5
#for D2
#vr=5
#trackvr=4
#
#for i in {1..310}
#for i in {90..240} #H2
#for i in {253..310} #D2
#for i in {90..140} #H2
#for i in {141..240} #H2
#for i in {176..176} #H2
#do
#for H2
#conf=conf/Run62/analyzer.hashi.20150723
#for D2
#conf=conf/Run62/calib.conf
#conf=conf/Run62/analyzer_run62_20150610.conf
#conf=/group/had/knucl/e15/kinoue/k18ana/conf/Run62/calib.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run62/myall/run62_%04d_v%02d.root $i $vr)
#in=$(printf data/Run62/run62_%04d.dat $i)
#track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
#done

#Run64
vr=5
trackvr=5
for i in {1..80}
do
conf=conf/Run64/analyzer_run64.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run64/myall/run64_%04d_v%02d.root $i $vr)
in=$(printf data/Run64/run64_%04d.dat $i)
track=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d.root $i $trackvr)
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
