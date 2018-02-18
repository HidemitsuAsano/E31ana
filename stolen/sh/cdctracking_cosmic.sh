#! /bin/sh

vr=12

#for i in {1..310} 
#for i in {115..151}  #CDC H.V./Vth study
for i in {115..163}
#for i in {176..180}
#for i in {253..265}
#for i in {266..270}
#for i in {271..275}
#for i in {276..285}
#for i in {286..290}
#for i in {291..300}
#for i in {209..209}
#for i in {253..310} #D2
#for i in {90..140} #H2
#for i in {141..240} #H2
do
#conf=conf/Run62/analyzer_run62_20150610.conf
conf=conf/Run62/analyzer_run62_cosmic.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run62/cdctracking/run62_%04d_v%02d_cosmic.root $i $vr)
in=$(printf data/Run62/run62_%04d.dat $i)
bsub -qs -o tmp.log ./bin/cdctracking_cosmic "$conf" "$evnum" "$out" "$in" 
done

