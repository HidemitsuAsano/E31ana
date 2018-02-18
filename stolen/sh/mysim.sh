#! /bin/sh

#Run65
rawvr=1
vr=1
#conf=conf/Run65/analyzer_run65.conf
conf=conf/Run65/analyzer_sim.conf
#conf=conf/Run65/tmp/analyzer_0005.conf
#outd=3NASpnpipipi
#outd=2NASpnpipi
#outd=2NALpn
#outd=3NALpnpi
#outd=3NASpn
#outd=3NApippn
outd=3NAS1385pn
#outd=2StepKdLp
#outd=2StepKppLp
#outd=2StepKpnLp
#outd=2StepS1385pLp
#outd=2StepL1405pLp
#outd=2StepL1520pLp
#outd=2StepLpLp
#outd=2StepSpLp
#outd=KppSp_40_50
#outd=3NASmpp
#outd=3NALpppim
#outd=3NALpppimpi
#outd=3NASpppim
#outd=3NASpppimpi

#outd=SLConv_S0
#outd=SLConv_Sp
#outd=SLConv_Sm

#out=$(printf root/mc/"$outd"/sim_%04d_v%02d.root $i $vr)
#in=$(printf /group/had/knucl/e15/kinoue/mc/k3He_Lpn/1/test_mc.root)
#bsub -ql -o tmp.log ./bin/mysim "$conf" "$out" "$in"

#for i in {1..1}
#for i in {1..100}
#for i in {1..500}
#for i in {1..1000}
for i in {1..1500}
#for i in {1..2000}
#for i in {501..1000}
#for i in {1001..1500}
do
out=$(printf tape/mc/"$outd"/sim_%04d_v%02d.root $i $vr)
in=$(printf tape/mcraw/"$outd"/sim_%04d_v%02d.root $i $rawvr)
bsub -ql -o tmp.log ./bin/mysim "$conf" "$out" "$in"
done
