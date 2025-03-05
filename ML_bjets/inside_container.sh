#!/bin/bash

dataname=$1 # The data set to run on
nJets=$2 # The number of jets used for training/testing
minpt=$3 # min pT of the selected jets
maxpt=$4 # max pT of the selected jets

export LD_LIBRARY_PATH="/home/hadi/HadiTest/Python/Python3/lib:$LD_LIBRARY_PATH"
export PATH="/home/hadi/HadiTest/Python/Python3/bin:$PATH"
export PYTHON_EXECUTABLE=/home/hadi/HadiTest/Python/Python3/bin/python #PATH Python interpreter executable
export PYTHON_INCLUDE_DIRS=/home/hadi/HadiTest/Python/Python3/include/Python3.4m #PATH Directory where to find Python.h
export PYTHON_LIBRARIES=/home/hadi/HadiTest/Python/Python3/lib #PATH Full path to Python library
export PYTHON_INCLUDE_DIR=/home/hadi/HadiTest/Python/Python3/include/Python3.4m #PATH Directory where to find Python.h
export PYTHON_LIBRARY=/home/hadi/HadiTest/Python/Python3/lib #PATH Full path to Python library
alias python=/home/hadi/HadiTest/Python/Python3/bin/Python3

source /home/hadi/HadiTest/root/bin/thisroot.sh

set -x

python3 --version

# echo $dataname $nJets $minpt $maxpt

python3 btagging.py $dataname $nJets $minpt $maxpt
