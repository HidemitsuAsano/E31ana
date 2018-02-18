#! /bin/sh

vr=3
trackvr=3

for i in {19..28} #USWK Scan with Empty target
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {79..88} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {153..161} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {239..247} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {282..290} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {345..353} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {425..433} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {501..509} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {561..569} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {661..669} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {712..720} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

for i in {771..779} #USWK Scan with 3He
do
conf=conf/Run65/analyzer_run65.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
in=$(printf data/Run65/run65_%04d.dat $i)
#track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
done

###################################################################
###################################################################
###################################################################

#
#for i in {114..114} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {136..136} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {260..260} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {281..281} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {423..424} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {444..444} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {469..469} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {493..500} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {534..534} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {559..560} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {583..583} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {608..608} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {635..635} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {659..660} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {686..686} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {711..711} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {745..745} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {769..770} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#for i in {819..820} #D5 Scan
#do
#conf=conf/Run65/analyzer_run65.conf
#evnum=$(printf %04d $i)
#out=$(printf root/ana/Run65/mypreana/run65_%04d_v%02d.root $i $vr)
#in=$(printf data/Run65/run65_%04d.dat $i)
##track=$(printf root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#track=$(printf /group/had/knucl/e15/noumi/root/ana/Run65/cdctracking/run65_%04d_v%02d.root $i $trackvr)
#bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in"
##bsub -qs -o tmp.log ./bin/mypreana "$conf" "$evnum" "$out" "$in" "$track"
#done
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
