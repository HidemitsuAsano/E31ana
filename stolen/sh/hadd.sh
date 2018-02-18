#! /bin/sh

vr=4

for i in {31..33} #pi scattering
do
in1=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d-1.root $i $vr)
in2=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d-2.root $i $vr)
in3=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d-3.root $i $vr)
in4=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d-4.root $i $vr)
in5=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d-5.root $i $vr)
out=$(printf root/ana/Run64/cdctracking/run64_%04d_v%02d.root $i $vr)
hadd "$out" "$in1" "$in2" "$in3" "$in4" "$in5"
done

