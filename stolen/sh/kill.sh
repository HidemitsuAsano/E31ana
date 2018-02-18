#! /bin/sh

vr=1


for i in {31763000..31763500}
do
out=$(printf %04d $i)
bkill "$out"
done

