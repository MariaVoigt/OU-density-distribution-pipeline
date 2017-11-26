#!/bin/bash

source /etc/profile.d/000-modules.sh

module purge
module load jdk/8 scala/2.11 ammonite

JAVA_OPTS=-Xmx32G \
	 command time -v amm \
	 ~/orangutan_density_distribution/src/validation/quartiles-categorized-diff.scala \
	 /work-local/maria/ \
	 /work/voigtma/bootstrap_no_0.1/bootstrap_diff_resource_use_grid.csv \
	 /data/idiv_kuehl/maria_data/bootstrap/resource_use_grid.csv \
	 1999,2015
