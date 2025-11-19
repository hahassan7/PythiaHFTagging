#!/bin/bash

condor_submit Submit.sub dataname=$1 nJets=$2 minpt=$3 maxpt=$4 outputname=$5
