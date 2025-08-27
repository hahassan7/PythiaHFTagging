#!/bin/bash

eval `alienv -w /home/hadi/alicesw/sw --no-refresh printenv ROOT/latest,fastjet/latest,pythia/latest,ONNXRuntime/latest`

set -x

# root.exe -l -b -q "${currentDir}/bjet_analysis.cpp(${NEVENTS},${PTHATMIN},${PTHATMAX})"

${currentDir}/bjet_analysis -n ${NEVENTS} -pmin ${PTHATMIN} -pmax ${PTHATMAX} -ml /scratch/project_2003583/focal-full-sim/HadiTest/PythiaHFTagging/ML_bjets/Models/pTHardBins/model.onnx
# ${currentDir}/bjet_analysis -n ${NEVENTS} -pmin ${PTHATMIN} -pmax ${PTHATMAX}

# sbatch submit.sh pThat_5_20 650000 5 20
# sbatch submit.sh pThat_20_40 630000 20 40
# sbatch submit.sh pThat_40_70 630000 40 70
# sbatch submit.sh pThat_70_1000 630000 70 -1

# sbatch submit.sh Eval_pThat_5_20 1000000 5 20
# sbatch submit.sh Eval_pThat_20_40 1000000 20 40
# sbatch submit.sh Eval_pThat_40_70 1000000 40 70
# sbatch submit.sh Eval_pThat_70_1000 650000 70 -1
