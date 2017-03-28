#!/bin/bash

#$ -N boot_abundPred

#$ -S /bin/bash

#$ -l h_rt=36:00:00
#$ -l h_vmem=50G,highmem

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log

#$ -j y

#$ -l avx

#$ -m ea


module load R
module load git 2> /dev/null

printf "current git version: %s" $(git rev-parse HEAD)
[[ -n $(git diff-index --name-only HEAD) ]] && echo "-dirty"
printf "\n"

mkdir -p ${OUTPUT_PATH:=/work/$USER/$JOB_NAME-$JOB_ID}

YEAR=$SGE_TASK_ID

Rscript \
    ~/orangutan_density_distribution/src/prediction/abundance_prediction_bootstrap_per_year.R \
    -o $OUTPUT_PATH \
    --year-to-predict $YEAR \
    --worker=${NSLOTS:-1} \
    "$@"

