 #! /bin/sh
cd /afs/cern.ch/work/o/oozcelik/GitTest/CMSSW_9_4_4/src/v0407/bmm4uml

eval `scramv1 runtime -sh`

echo 'running the code'

for number in {1..4}
do
echo "runnig the set of toy =  $number"
./bmm4toystudy genfit minos tausplot seed=1214${number} iterations=50. yield_scale=200.0 ws_gen=wspace_prefit.root ws_res=wspace_toyresult_200_v4_${number}.root
done
