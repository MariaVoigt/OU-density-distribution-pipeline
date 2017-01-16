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

mkdir -p ${OUTPUT_PATH:=/work/$USER/$JOB_NAME-$JOB_ID}

YEAR=$SGE_TASK_ID

export MC_CORES=${NSLOTS:-1}

Rscript \
    ~/orangutan_density_distribution/src/prediction/abundance_prediction.R \
    -o $OUTPUT_PATH \
    --year-to-predict $YEAR \
    "$@"

