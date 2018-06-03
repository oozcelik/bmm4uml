#!/bin/bash
cd /afs/cern.ch/work/o/oozcelik/GitTest/CMSSW_9_4_4/src/bmm4uml
eval `scramv1 runtime -sh`

dcacheDirIn="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/ursl/bmm4/small/180428/s01/"
filestocopy=""


for files in `gfal-ls $dcacheDirIn | grep 2016`
do
  echo $files
  filestocopy=root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/ursl/bmm4/small/180428/s01/$files;
  xrdcp $filestocopy input/bmm4/
done
