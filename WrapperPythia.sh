#!/bin/bash

mkdir -p /afs/cern.ch/work/h/hahassan/PythiaHFTagging/OutputFiles/Default_settings/${2}_${3}
condor_submit SubmitPythia.sub NEVENTS=$1 PTHATMIN=$2 PTHATMAX=$3
