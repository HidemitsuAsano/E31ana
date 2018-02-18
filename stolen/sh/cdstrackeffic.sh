#! /bin/sh

#vr=12
#travr=12
#conf=conf/Run62/analyzer_run62_CDCHVScan.conf
#for i in {115..163}  #CDC H.V./Vth study
#for i in {115..151}  #CDC H.V./Vth study
#vr=12
#travr=5
#conf=conf/Run62/analyzer.hashi.20150723
#for i in {153..163}  #a part of H2 run

#for i in {26..43}
#for i in {44..55}
#for i in {1..310}
vr=7
travr=7
conf=conf/Run64/analyzer_run64.conf
for i in {31..33}  #CDC H.V./Vth study
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run64/cds/run64_%04d_v%02d_trackeffic.root $i $vr)
in=$(printf data/Run64/run64_%04d.dat $i)
track=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d.root $i $travr)
bsub -ql -o tmp.log ./bin/cdstrackeffic "$conf" "$evnum" "$out" "$in" "$track"
done
