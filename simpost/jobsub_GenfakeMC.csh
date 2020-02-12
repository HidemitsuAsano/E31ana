#!/bin/tcsh -f

bsub -q s root -q -b GenPiPinFake.C+\(1\)
bsub -q s root -q -b GenPiPinFake.C+\(2\)
bsub -q s root -q -b GenPiPinFake.C+\(3\)
bsub -q s root -q -b GenPiPinFake.C+\(4\)
bsub -q s root -q -b GenPiPinFake.C+\(5\)
bsub -q s root -q -b GenPiPinFake.C+\(6\)
bsub -q s root -q -b GenPiPinFake.C+\(7\)
bsub -q s root -q -b GenPiPinFake.C+\(8\)
bsub -q s root -q -b GenPiPinFake.C+\(9\)
bsub -q s root -q -b GenPiPinFake.C+\(10\)


set hname=`hostname -s`
echo ${hname}
while (1)
  @ njob=`bjobs | grep ${hname} | wc -l`
  if ( $njob < 1 ) then 
    hadd -f fakepippim_pippimn_sum.root fakepippim_pippimn[1-9].root fakepippim_pippimn10.root
    break
  endif
  sleep 60
end
