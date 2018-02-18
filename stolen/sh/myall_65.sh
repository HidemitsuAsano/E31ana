#! /bin/sh

#Run65
vr=1
trackvr=7
conf=conf/Run65/analyzer_run65.conf
confe=conf/Run65/analyzer_run65_Empty.conf
#conf=conf/Run65/tmp/analyzer_run65_0002.conf
#confe=conf/Run65/analyzer_run65_Empty.conf

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

for i in {30..41} #Empty
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$confe" "$evnum" "$out" "$in" "$track"
done
for i in {101..152} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {164..208} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {248..280} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {291..344} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {358..422} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {435..492} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {510..558} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {571..658} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {671..710} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {721..768} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
for i in {780..846} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in"
bsub -ql -o tmp.log ./bin/myall "$conf" "$evnum" "$out" "$in" "$track"
done
