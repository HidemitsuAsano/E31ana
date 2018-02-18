#! /bin/sh

#Run65
vr=1
trackvr=7
#trackvr=1
anavr=1
conf=conf/Run65/analyzer_run65.conf
confe=conf/Run65/analyzer_run65_Empty.conf

#for i in {30..41} #Empty
#do
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
##track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -ql -o tmp.log ./bin/myreall "$confe" "$evnum" "$out" "$in" "$track" "$all"
#done

for i in {101..152} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {164..208} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {248..280} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {291..344} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {358..422} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {435..492} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {510..558} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {571..658} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {671..710} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {721..768} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

for i in {780..846} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/myreall/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
all=$(printf root/ana/Run65/myall/run65_%04d_v%02d.root $i $anavr)
bsub -ql -o tmp.log ./bin/myreall "$conf" "$evnum" "$out" "$in" "$track" "$all"
done

