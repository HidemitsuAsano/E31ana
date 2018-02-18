#! /bin/sh

vr=6

#Run62
for i in {78..89} #Empty
#for i in {121..160} #3He
#for i in {101..324} #3He
#for i in {101..399} #3He
do
conf=conf/Run62/analyzer_run62.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
done

