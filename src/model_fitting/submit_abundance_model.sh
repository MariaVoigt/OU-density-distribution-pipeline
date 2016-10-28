#!/bin/bash

#$ -N abundMod

#$ -S /bin/bash

#$ -l h_rt=720:00:00
#$ -l h_vmem=18G,highmem

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -pe smp 20


#$ -m ea

module load R

if [[ -z $1 ]] ; then
  echo "qsub $0 /path/to/input"
  exit 1
fi

if [[ ! -d $1 ]] ; then
  echo "directory does not exist: $1"
  exit 1
fi

INPUT_PATH=$1
OUTPUT_PATH=/work/$USER/$JOB_NAME-$JOB_ID
FUN_FILE=$2


mkdir $OUTPUT_PATH

export MC_CORES=${NSLOTS:-1}

Rscript \
$HOME/model_and_predictions/abundance_model.R \
    $INPUT_PATH \
    $OUTPUT_PATH \
    $FUN_FILE

