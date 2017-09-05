#!/bin/bash

# ------------------------------------------------------------------------------
# config
# ------------------------------------------------------------------------------


# here somehow write the prefix depending on what we are testing
JOB_NAME_PREFIX=ppln_ae75m_50

INPUT_PATH='/data/idiv_kuehl/maria_data/pipeline'
NAME=$JOB_NAME_PREFIX-$(date +%FT%H-%M-%S)

# FIX THIS
OUTPUT_PATH=/work/$USER/$NAME
# qsub -v means that OUTPUT_PATH goes into environment (normally export, but with qsub different)
QSUB="qsub -terse -v OUTPUT_PATH=$OUTPUT_PATH -o /work/$USER/$NAME.log -j y"


# ------------------------------------------------------------------------------
# model-fitting & prediction 
# ------------------------------------------------------------------------------

MODEL_JOB_ID_FITTING=$($QSUB \
			   -N ${JOB_NAME_PREFIX}_fitting \
			   $HOME/orangutan_density_distribution/src/model_fitting/submit_abundance_model.sh \
			   -i $INPUT_PATH \
			   --ESW-aerial 0.075 \
			   --include-aerial \
			   --stability)


  # prediction submitten
 MODEL_JOB_ID_PRED=$($QSUB \
			 -N ${JOB_NAME_PREFIX}_prediction \
			 -hold_jid $MODEL_JOB_ID_FITTING \
			 -t 1999:2015 \
			 $HOME/orangutan_density_distribution/src/prediction/submit_abundance_pred.sh \
			 -i $INPUT_PATH \
			 --input-model $OUTPUT_PATH \
			 --focal-change-predictor NA)

