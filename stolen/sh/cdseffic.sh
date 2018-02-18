#! /bin/sh

vr=7
#layer=1

#for i in {26..43}
#for i in {44..55}
#for i in {1..310}
#for i in {115..151}  #CDC H.V./Vth study
#do
#for layer in {2..15}
#do
##conf=conf/MyStudy/analyzer_run62.conf
##conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/Run62/analyzer_run62_CDCHVScan.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run62/cds/run62_%04d_v%02d_effic%02d.root $i $vr $layer)
#in=$(printf data/Run62/run62_%04d.dat $i)
#track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d_kill%02d.root $i $vr $layer)
#bsub -qs -o tmp.log ./bin/cdseffic "$conf" "$evnum" "$out" "$in" "$track"
#done
#done
for i in {31..33}  #CDC H.V./Vth study
do
for layer in {1..15}
do
#conf=conf/MyStudy/analyzer_run62.conf
#conf=conf/Run62/analyzer_run62_20150610.conf
conf=conf/Run64/analyzer_run64.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run64/cds/run64_%04d_v%02d_effic%02d.root $i $vr $layer)
in=$(printf data/Run64/run64_%04d.dat $i)
track=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d_kill%02d.root $i $vr $layer)
bsub -qs -o tmp.log ./bin/cdseffic "$conf" "$evnum" "$out" "$in" "$track"
done
done
