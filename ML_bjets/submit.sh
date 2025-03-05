#!/bin/bash
#SBATCH --job-name=bjetML
#SBATCH --account=project_2003583
##SBATCH --begin=now+3hour
#SBATCH --partition=gpu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16000
#SBATCH --gres=gpu:v100:1,nvme:10

#SBATCH --array=1-1 #defines SLURM_ARRAY_TASK_ID
# If you change output/error here, please change
# the mv command at the end of this macro
#SBATCH --output=logs/output_%A-run%a.txt
#SBATCH --error=logs/errors_%A-run%a.txt

dataname=$1 # The data set to run on
nJets=$2 # The number of jets used for training/testing
minpt=$3 # min pT of the selected jets
maxpt=$4 # max pT of the selected jets


n=$SLURM_ARRAY_TASK_ID

echo "Starting Slurm array job ${SLURM_ARRAY_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}"

/scratch/project_2003583/focal-full-sim/HadiTest/ML_bjets/run $dataname $nJets $minpt $maxpt
sleep 1

# --mem-per-gpu=16000
