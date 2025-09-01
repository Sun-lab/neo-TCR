#!/bin/bash

# M13
cd /fh/fast/sun_w/yub/tcr/tracerx/tcr-seq/decombinator_results_M13
sbatch \
  --job-name dcr_M13 \
  --cpus-per-task 1 \
  --mem 4G \
  --array=1-$(wc -l < Batch_M13.txt) \
  --wrap='cmd=$(sed -n ${SLURM_ARRAY_TASK_ID}p Batch_M13.txt);
          eval "$cmd"'

# I8
cd /fh/fast/sun_w/yub/tcr/tracerx/tcr-seq/decombinator_results_I8
sbatch \
  --job-name dcr_I8 \
  --cpus-per-task 1 \
  --mem 4G \
  --array=1-$(wc -l < Batch_I8.txt) \
  --wrap='cmd=$(sed -n ${SLURM_ARRAY_TASK_ID}p Batch_I8.txt);
          eval "$cmd"'