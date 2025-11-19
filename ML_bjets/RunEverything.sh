#!/bin/sh
date

dataname=$1 # The data set to run on
nJets=$2 # The number of jets used for training/testing
minpt=$3 # min pT of the selected jets
maxpt=$4 # max pT of the selected jets
outputname=$5 # The name of the output folder

echo "Executing run $7 on" `hostname` "on $dataname with $nJets jets and pT range [$minpt,$maxpt]"

cd /afs/cern.ch/work/h/hahassan/PythiaHFTagging/ML_bjets/

source /afs/cern.ch/work/h/hahassan/env.sh

set -x

python3 --version

# echo $dataname $nJets $minpt $maxpt

python3 -c "import tensorflow as tf; print(\"GPU available\" if tf.config.list_physical_devices('GPU') else \"GPU not available\")"

time python3 btagging.py $dataname $nJets $minpt $maxpt $outputname
