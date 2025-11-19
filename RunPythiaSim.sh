#!/bin/bash

export NEVENTS=${3}
export PTHATMIN=${4}
export PTHATMAX=${5}

eval $(/cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@O2Physics::daily-20251113-0000-1)

set -x

time ./bjet_analysis -n ${NEVENTS} -pmin ${PTHATMIN} -pmax ${PTHATMAX} -ml /afs/cern.ch/work/h/hahassan/PythiaHFTagging/ML_bjets/Models/pTHardBins_TrainValTest_Split/model.onnx
