#! /bin/sh

vr=7
conf=conf/Run65/analyzer_run65.conf

#Run65
#for i in {29..41} #Empty
#for i in {101..101} #3He
#for i in {101..208} #3He
#for i in {209..280} #3He
#for i in {281..344} #3He
#for i in {345..422} #3He
#for i in {423..469} #3He
#for i in {470..491} #3He
#for i in {1..846} #3He

#for i in {101..301} #3He
#for i in {302..551} #3He
#for i in {552..846} #3He

#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
#bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
##bsub -ql -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#done

#2016.06.29 DONE ------------>
#for i in {164..208} #3He
#do
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
#bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#done
#for i in {248..280} #3He
#do
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
#bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#done
#<--------- 2016.06.29 DONE

#2016.07.01 DONE ------------>
#for i in {435..492} #3He
#do
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $vr)
#in=$(printf /group/had/knucl/e15/data/Run65/run65_%04d.dat $i)
#bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
##bsub -ql -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#done
#for i in {510..558} #3He
#do
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $vr)
#in=$(printf /group/had/knucl/e15/data/Run65/run65_%04d.dat $i)
#bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
##bsub -ql -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#done
#<--------- 2016.07.01 DONE

#2016.07.03 DONE ------------>
for i in {721..768} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $vr)
in=$(printf /group/had/knucl/e15/data/Run65/run65_%04d.dat $i)
bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#bsub -ql -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
done
for i in {780..846} #3He
do
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $vr)
in=$(printf /group/had/knucl/e15/data/Run65/run65_%04d.dat $i)
bsub -qh -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
#bsub -ql -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
done
#<--------- 2016.07.03 DONE
