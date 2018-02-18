#! /bin/sh

vr=7

#for i in {1..310} 
#for i in {115..151}  #CDC H.V./Vth study
#for i in {115..163}  #CDC H.V./Vth study
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
#Run64
#for i in {1..80} #pi scattering
#for i in {31..33} #pi scattering
for i in {81..105} #pi scattering
do
#conf=conf/Run62/analyzer_run62_20150610.conf
#conf=conf/Run62/analyzer_run62_CDCHVScan.conf
conf=conf/Run64/analyzer_run64.conf
evnum=$(printf %04d $i)
out=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d.root $i $vr)
in=$(printf data/Run64/run64_%04d.dat $i)
bsub -ql -o tmp.log ./bin/cdctracking "$conf" "$evnum" "$out" "$in" 
done

