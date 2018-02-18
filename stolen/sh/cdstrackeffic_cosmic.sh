#! /bin/sh

vr=12

#for i in {26..43}
#for i in {44..55}
#for i in {1..310}
for i in {115..163}  #CDC H.V./Vth study
do
#conf=conf/MyStudy/analyzer_run62.conf
#conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/Run62/analyzer_run62_CDCHVScan.conf
conf=conf/Run62/analyzer_run62_cosmic.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/cds/run62_%04d_v%02d_trackeffic_cosmic.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d_cosmic.root $i $vr)
bsub -qs -o tmp.log ./bin/cdstrackeffic_cosmic "$conf" "$evnum" "$out" "$in" "$track"
done
