#! /bin/sh

vr=1

for i in {242..250}
do
conf=conf/MyStudy/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/bhd/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
nohup ./bin/bhd "$conf" "$evnum" "$out" "$in" &
done

#renice -n 15 -u yamaga

