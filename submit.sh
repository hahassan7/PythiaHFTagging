#!/bin/bash
## sbatch usage: sbatch submit.sh <comment> <generator> <nevents> <geometry>
## sbatch will check arguments from the comments in the
## beginning of this file.
#SBATCH --job-name=PYTHIA8
#SBATCH --account=project_2003583
#SBATCH --partition=small
##SBATCH --begin=now+3hour
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-20 #defines SLURM_ARRAY_TASK_ID
# If you change output/error here, please change
# the mv command at the end of this macro
#SBATCH --output=logs/output_%A-run%a.txt
#SBATCH --error=logs/errors_%A-run%a.txt

comment=$1 # Identifying different runs
nevents=$2
pTHatMin=$3 # Min pT hat for the hard scattering
pTHatMax=$4 # Max pT hat for the hard scattering

n=$SLURM_ARRAY_TASK_ID
outputdir=/scratch/project_2003583/focal-full-sim/HadiTest/PythiaHFTagging/GeneratedEvents/${comment}/$(printf "%03d" $n)

echo "Starting Slurm array job ${SLURM_ARRAY_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}"

mkdir -p $outputdir
mkdir -p ${outputdir}/logs
/scratch/project_2003583/focal-full-sim/HadiTest/PythiaHFTagging/run $outputdir $nevents $pTHatMin $pTHatMax
sleep 1
mv logs/output_${SLURM_ARRAY_JOB_ID}-run${SLURM_ARRAY_TASK_ID}.txt ${outputdir}/logs/
mv logs/errors_${SLURM_ARRAY_JOB_ID}-run${SLURM_ARRAY_TASK_ID}.txt ${outputdir}/logs/
sleep 1

