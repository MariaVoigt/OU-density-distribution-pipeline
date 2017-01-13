#!/bin/bash

#$ -N abundMod

#$ -S /bin/bash

#$ -l h_rt=48:00:00
#$ -l h_vmem=3G

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -pe smp 10-28


#$ -m ea

#$ -cwd

module load R
module load git

printf "current git version: %s" $(git rev-parse HEAD)
[[ -n $(git diff-index --name-only HEAD) ]] && echo "-dirty"
printf "\n"

OUTPUT_PATH=/work/$USER/$JOB_NAME-$JOB_ID

mkdir -p $OUTPUT_PATH

export MC_CORES=${NSLOTS:-1}

Rscript \
    abundance_model.R \
    -o $OUTPUT_PATH \
    "$@"