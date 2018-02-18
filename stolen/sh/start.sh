#! /bin/sh

vr=1

#for i in {31..41}
#for i in {42..50}
#for i in {51..60}
#for i in {61..70}
#for i in {71..82}
for i in {201..210}
do
conf=conf/MyStudy/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/start/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
nohup ./bin/start "$conf" "$evnum" "$out" "$in" &
done

renice -n 15 -u yamaga

