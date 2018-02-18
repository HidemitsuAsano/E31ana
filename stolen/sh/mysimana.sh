#! /bin/sh

#Run65
rawvr=1
simvr=1
vr=3
conf=conf/Run65/analyzer_sim.conf
#outd=2NALpn
#outd=2NALpnpi
#outd=2NALpnpipi
#outd=2NALpnpipipi
#outd=2NASpn
#outd=2NASpnpi
#outd=2NASpnpipi
#outd=2NASpnpipipi

outd=3NALpn
#outd=3NALpnpi
#outd=3NALpnpipi
#outd=3NALpnpipipi
#outd=3NASpn
#outd=3NASpnpi
#outd=3NASpnpipi
#utd=3NASpnpipipi

#outd=3NALpppim
#outd=3NALpppimpi

#outd=3NASpppim
#outd=3NASpppimpi

#outd=3NAS1385pn
#outd=3NAS1385pnpi

#outd=3NApippn
#outd=3NASmpp

#outd=2StepKdLp
#outd=2StepKppLp
#outd=2StepKpnLp
#outd=2StepS1385pLp
#outd=2StepL1405pLp
#outd=2StepL1520pLp
#outd=2StepLpLp
#outd=2StepSpLp

#outd=SLConv_S0
#outd=SLConv_Sp
#outd=SLConv_Sm


#outd=KppLp_40_50

#out=$(printf root/mc/k3He_S0pn/sim_%04d_v%02d.root $i $vr)
#in=$(printf /group/had/knucl/e15/kinoue/mc/k3He_S0pn/%d/test_mc.root $i)

#for i in {1..20}
#for i in {1..100}
for i in {1..500}
#for i in {1..1000}
#for i in {1..1500}
#for i in {1..2000}
#for i in {501..1000}
do
out=$(printf root/mcana/"$outd"/sim_%04d_v%02d.root $i $vr)
raw=$(printf tape/mcraw/"$outd"/sim_%04d_v%02d.root $i $rawvr)
sim=$(printf tape/mc/"$outd"/sim_%04d_v%02d.root $i $simvr)
bsub -ql -o tmp.log ./bin/mysimana "$conf" "$out" "$raw" "$sim"
done
