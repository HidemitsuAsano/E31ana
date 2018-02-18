#! /bin/sh

vr=12
#vr=5

#for i in {26..43}
#for i in {44..55}
#for i in {1..310}
#for i in {115..151}  #CDC H.V./Vth study
#for i in {153..162}  #CDC H.V./Vth study
for i in {115..162}  #CDC H.V./Vth study
do
#conf=conf/MyStudy/analyzer_run62.conf
#conf=conf/Run62/analyzer_run62_20150610.conf
conf=conf/Run62/analyzer_run62_CDCHVScan.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/cds/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $vr)
bsub -qs -o tmp.log ./bin/cds "$conf" "$evnum" "$out" "$in" "$track"
done
