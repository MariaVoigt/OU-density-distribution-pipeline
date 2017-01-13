#!/bin/bash

#$ -N abundPred

#$ -S /bin/bash

#$ -l h_rt=48:00:00
#$ -l h_vmem=3G

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log

#$ -j y

#$ -pe smp 10-28
#$ -l avx

#$ -m ea


module load R
module load git

printf "current git version: %s" $(git rev-parse HEAD)
[[ -n $(git diff-index --name-only HEAD) ]] && echo "-dirty"
printf "\n"

if [[ ! -d $1 ]] ; then
  echo "directory does not exist: $1"
  exit 1
fi

INPUT_PATH=$1
OUTPUT_PATH=/work/$USER/$JOB_NAME-$JOB_ID
FUN_FILE=$2
PRED_FILE=$3
YEAR=$SGE_TASK_ID

mkdir -p $OUTPUT_PATH

export MC_CORES=${NSLOTS:-1}

Rscript \
    ~/orangutan_density_distribution/src/prediction/abundance_prediction.R \
    $INPUT_PATH \
    $OUTPUT_PATH \
    $FUN_FILE \
    $PRED_FILE \
    $YEAR

