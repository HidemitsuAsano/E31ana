#! /bin/sh

vr=1

#for i in {9..20}
for i in {20..33}
do
#if [ $i -ne 10 ] ; then
conf=conf/MyStudy/fdc_analyzer.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run61pre/fdc/FDCrun%04d_v%02d.root $i $vr)
in=$(printf data/Run61pre/FDC/FDCrun%04d.dat $i)
nohup ./bin/fdc "$conf" "$evnum" "$out" "$in" &
#fi
done

renice -n 15 -u yamaga

