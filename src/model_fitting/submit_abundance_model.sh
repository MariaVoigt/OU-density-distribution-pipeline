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

shift

mkdir $OUTPUT_PATH

export MC_CORES=${NSLOTS:-1}

Rscript \
    abundance_model.R \
    -i $INPUT_PATH \
    -o $OUTPUT_PATH \
    "$@"
