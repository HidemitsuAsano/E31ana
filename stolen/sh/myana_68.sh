#! /bin/sh

#Run68
vr=1
trackvr=2
anavr=1
conf=conf/Run68/analyzer_run68.conf
for i in {49..90}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {102..112}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {125..140}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {152..157}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {169..177}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {190..200}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {212..236}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {248..272}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {292..314} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {319..334} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {337..352} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
for i in {355..359} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/myana/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
all=$(printf root/ana/Run68/myall/run68_%04d_v%02d.root $i $anavr)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/myana "$conf" "$evnum" "$out" "$in" "$track" "$all"
done
