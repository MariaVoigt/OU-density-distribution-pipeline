#!/bin/bash

source /etc/profile.d/000-modules.sh

module purge
module load jdk/8 scala/2.11 ammonite

JAVA_OPTS=-Xmx32G \
parallel \
  --line-buffer \
  -j10 \
  'command time -v amm ~/orangutan_density_distribution/src/validation/quartiles-categorized.scala /work-local/maria/ /work/voigtma/bootstrap_no_0.1_resource_use/bootstrap_abundance_{/.}_output.csv {} 1999,2015' ::: /data/idiv_kuehl/maria_data/bootstrap/resource_use_grid.csv
