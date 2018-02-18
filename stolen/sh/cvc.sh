#! /bin/sh

vr=1

for i in {48..48}
do
#if [ $i -ne 10 ] ; then
conf=conf/MyStudy/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/cvc/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
nohup ./bin/cvc "$conf" "$evnum" "$out" "$in" &
#fi
done

renice -n 15 -u yamaga

