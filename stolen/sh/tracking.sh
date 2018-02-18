#! /bin/sh

vr=9

for i in {1..1}
do
#if [ $i -ne 3 ] ; then
conf=conf/MyStudy/cdc_analyzer.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run61pre/tracking/CDCrun%04d_v%02d.root $i $vr)
in=$(printf data/Run61pre/CDC/CDCrun%04d.dat $i)
nohup ./bin/evtracking "$conf" "$evnum" "$out" "$in" &
#echo "./bin/evtracking "$conf" "$evnum" "$out" "$in""
#fi
done

renice -n 15 -u yamaga

