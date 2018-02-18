#! /bin/sh

#Run65
vr=90
trackvr=7
conf=conf/Run65/analyzer_run65.conf
confe=conf/Run65/analyzer_run65_Empty.conf

#for i in {101..111}
#for i in {101..500}
#for i in {501..846}
#for i in {101..160}
#for i in {161..260}
#for i in {101..846}
#for i in {79..88}
#for i in {101..469} #3He (Draw = +983A)
#for i in {209..469} #3He (Draw = +983A)
#for i in {470..820} #3He (Draw = -983A)
#for i in {821..846} #3He (Draw = +983A)
#for i in {470..491} #3He (Draw = -983A)
#for i in {470..846}
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
##track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qh -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
#done

for i in {185..185}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {583..583}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {635..635}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {534..534}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {136..136}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {444..444}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {819..819}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {260..260}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {343..343}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {686..686}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {608..608}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {469..469}
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
