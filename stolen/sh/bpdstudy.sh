#! /bin/sh

vr=1

for i in {42..42}
do
if [ $i -ne 10 ] ; then
conf=conf/MyStudy/bpd_analyzer.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run61pre/bpdstudy/BPDStudy%04d_v%02d.root $i $vr)
in=$(printf data/Run61pre/BPD/BPDrun%04d.dat $i)
nohup ./bin/bpdstudy "$conf" "$evnum" "$out" "$in" &
#echo "./bin/evtracking "$conf" "$evnum" "$out" "$in""
fi
done

#renice -n 15 -u yamaga

