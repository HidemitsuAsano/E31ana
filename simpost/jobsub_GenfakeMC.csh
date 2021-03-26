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
bsub -q s root -q -b GenPiPinFake.C+\(11\)
bsub -q s root -q -b GenPiPinFake.C+\(12\)
bsub -q s root -q -b GenPiPinFake.C+\(13\)
bsub -q s root -q -b GenPiPinFake.C+\(14\)
bsub -q s root -q -b GenPiPinFake.C+\(15\)
bsub -q s root -q -b GenPiPinFake.C+\(16\)
bsub -q s root -q -b GenPiPinFake.C+\(17\)
bsub -q s root -q -b GenPiPinFake.C+\(18\)
bsub -q s root -q -b GenPiPinFake.C+\(19\)
bsub -q s root -q -b GenPiPinFake.C+\(20\)


set hname=`hostname -s`
echo ${hname}
while (1)
  @ njob=`bjobs | grep ${hname} | wc -l`
  if ( $njob < 1 ) then 
    cd /gpfs/group/had/knucl/e15/asano/sim/fakemc_test/
    hadd -f fakepippim_pippimn_sum.root fakepippim_pippimn[1-9].root fakepippim_pippimn1[0-9].root fakepippim_pippimn20.root
    break
  endif
  sleep 60
end
