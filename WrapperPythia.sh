#!/bin/bash

mkdir -p /afs/cern.ch/work/h/hahassan/PythiaHFTagging/OutputFiles/Testing/${2}_${3}
# mkdir -p /afs/cern.ch/work/h/hahassan/PythiaHFTagging/OutputFiles/TreeAnalysisResults/${2}_${3}
# mkdir -p /afs/cern.ch/work/h/hahassan/PythiaHFTagging/OutputFiles/HFTaggingTree/${2}_${3}
condor_submit SubmitPythia.sub NEVENTS=$1 PTHATMIN=$2 PTHATMAX=$3
