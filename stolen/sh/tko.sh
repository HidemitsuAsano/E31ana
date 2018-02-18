#! /bin/sh

vr=1

#for i in {31..41}
#for i in {42..50}
#for i in {51..60}
#for i in {61..70}
#for i in {71..82}
for i in {51..100}
do
#if [ $i -ne 10 ] ; then
conf=conf/MyStudy/analyzer_timeshift.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62pre/timeshift/TimeShift%04d_v%02d.root $i $vr)
in=$(printf data/Run62pre/TimeShift%04d.dat $i)
nohup ./bin/start "$conf" "$evnum" "$out" "$in" &
#conf=conf/MyStudy/analyzer_run62.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run62pre/tko/TDCCalib%04d_v%02d.root $i $vr)
#in=$(printf data/Run62pre/TDCCalib/TDCCalib%04d.dat $i)
#nohup ./bin/evanatko "$conf" "$evnum" "$out" "$in" &
#fi
done

renice -n 15 -u yamaga

