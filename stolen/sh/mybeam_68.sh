#! /bin/sh

vr=1
trackvr=2
conf=conf/Run68/analyzer_run68.conf

for i in {49..90} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {102..112} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {125..140} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {152..157} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {159..166} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {169..177} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {190..200} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {212..236} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {248..272} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {292..314} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {319..334} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {337..352} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {355..359} #D2 Production
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/beam/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
#track=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/kinoue/Run68/evtracking_%d.root $i)
bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/beam "$conf" "$evnum" "$out" "$in" "$track"
done
