#!/bin/bash

eval `alienv -w /home/hadi/alicesw/sw --no-refresh printenv ROOT/latest,fastjet/latest,pythia/latest`

set -x

root.exe -l -b -q "${currentDir}/bjet_analysis.cpp(${NEVENTS},${PTHATMIN},${PTHATMAX})"

# ${currentDir}/bjet_analysis -n ${NEVENTS} -pmin ${PTHATMIN} -pmax ${PTHATMAX}
