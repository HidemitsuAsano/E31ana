for i in {41..44}
do 
echo $i
./bin/evanaraw conf/Run74/analyzer.conf $i run74raw_00$i.root /w/e15/data/Run74/run74_00$i.dat 
done
