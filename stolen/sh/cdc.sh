#! /bin/sh

vr=12

#for i in {1..310}
#for i in {115..151}  #CDC H.V./Vth study
for i in {153..163}
do
#conf=conf/Run62/analyzer_run62_20150610.conf
conf=conf/Run62/analyzer_run62_CDCHVScan.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/cdc/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
bsub -qs -o tmp.log ./bin/cdc "$conf" "$evnum" "$out" "$in" 
done

