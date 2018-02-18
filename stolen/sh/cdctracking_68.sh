#! /bin/sh

vr=2

#Run68
#for i in {49..65} #D2 Production
#for i in {71..82} #D2 Production
#for i in {83..91} #D2 Production
#for i in {102..112} #D2 Production
#for i in {161..162} #D2 Production
#for i in {177..177} #D2 Production
#for i in {190..193} #D2 Production
#for i in {194..195} #D2 Production
for i in {212..230} #D2 Production

do
conf=conf/Run68/analyzer_run68.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run68/cdctracking/run68_%04d_v%02d.root $i $vr)
in=$(printf data/Run68/run68_%04d.dat $i)
bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#bsub -ql -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
done

