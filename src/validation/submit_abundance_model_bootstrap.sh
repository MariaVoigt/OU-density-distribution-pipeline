#!/bin/bash

#$ -N boot

#$ -S /bin/bash

#$ -l h_rt=168:00:00
#$ -l h_vmem=6G

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -pe smp 20-28


#$ -m ea


module load R
module load git 2> /dev/null

printf "current git version: %s" $(git rev-parse HEAD)
[[ -n $(git diff-index --name-only HEAD) ]] && echo "-dirty"
printf "\n"

mkdir -p ${OUTPUT_PATH:=/work/$USER/$JOB_NAME-$JOB_ID}

export MC_CORES=${NSLOTS:-1}

Rscript \
    ~/orangutan_density_distribution/src/validation/marias_boot_correct.R \
    -o $OUTPUT_PATH \
    "$@"
