#! /bin/sh

#vr=12
#travr=12
#conf=conf/Run62/analyzer_run62_CDCHVScan.conf
#for i in {115..151}  #CDC H.V./Vth study
vr=12
travr=5
conf=conf/Run49c/analyzer.cosmic
for i in {100..110}  #a part of H2 run

#for i in {26..43}
#for i in {44..55}
#for i in {1..310}
do
#conf=conf/MyStudy/analyzer_run62.conf
#conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/Run62/analyzer_run62_CDCHVScan.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run49c/cds/run49c_%04d_v%02d_trackeffic_cosmic.root $i $vr)
in=$(printf data/Run49c/run49c_%04d.dat $i)
track=$(printf /home/had/thashi/grp/e15/hashimoto/run49c/cosmic_20130822_%d.root $i)
bsub -ql -o tmp.log ./bin/cdstrackeffic_cosmic "$conf" "$evnum" "$out" "$in" "$track"
done
