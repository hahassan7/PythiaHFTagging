#!/bin/bash

eval `alienv -w /home/hadi/alicesw/sw --no-refresh printenv ROOT/latest,fastjet/latest,pythia/latest,ONNXRuntime/latest`

set -x

# root.exe -l -b -q "${currentDir}/bjet_analysis.cpp(${NEVENTS},${PTHATMIN},${PTHATMAX})"

# ${currentDir}/bjet_analysis -n ${NEVENTS} -pmin ${PTHATMIN} -pmax ${PTHATMAX} -ml /scratch/project_2003583/focal-full-sim/HadiTest/PythiaHFTagging/ML_bjets/Models/pTHardBins/model.onnx
${currentDir}/bjet_analysis -n ${NEVENTS} -pmin ${PTHATMIN} -pmax ${PTHATMAX}
