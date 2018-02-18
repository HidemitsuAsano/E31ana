#! /bin/sh

vr=1
trackvr=4
anavr=1

for i in {102..114}
do
conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/MyStudy/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/myana2/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
ana=$(printf root/ana/Run62/myana/run62_%04d_v%02d.root $i $anavr)
bsub -qe -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$ana"
done
for i in {176..208}
do
conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/MyStudy/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/myana2/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
ana=$(printf root/ana/Run62/myana/run62_%04d_v%02d.root $i $anavr)
bsub -qe -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$ana"
done

for i in {211..240}
do
conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/MyStudy/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/myana2/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
track=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $trackvr)
ana=$(printf root/ana/Run62/myana/run62_%04d_v%02d.root $i $anavr)
bsub -qe -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$ana"
done

