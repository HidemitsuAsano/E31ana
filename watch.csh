#!/bin/tcsh -f
set hname=`hostname -s`
echo $hname
while (1)
#  @ njob=`bjobs | grep ${hname} | wc -l`
  @ njob=`bjobs | wc -l`
  @ nrunjob = `bjobs | grep RUN | wc -l`
  @ nrunjobs = `bjobs | grep "RUN   s" | wc -l`
  @ nrunjobl = `bjobs | grep "RUN   l" | wc -l`
  @ npendjobs = `bjobs | grep "PEND  s" | wc -l`
  @ npendjobl = `bjobs | grep "PEND  l" | wc -l`
  set prit = `bqueues -l s | grep hiasano`
  echo "job priority $prit[3]"
#  bqueues -l s | grep hiasano | awk '{print $3}'
  echo "$nrunjob/$njob jobs running" 
  echo "sjob $nrunjobs, pend $npendjobs" 
  echo "ljob $nrunjobl, pend $npendjobl" 
  gquota -G /group/had/knucl
  
  sleep 60
end
