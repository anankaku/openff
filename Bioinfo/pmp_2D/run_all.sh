#!/bin/bash

mkdir -p slurm_logs

shopt -s nullglob
for file in *.com; do
  base="${file%.com}"
  echo "Submitting $file"
  sbatch --job-name="$base" run_single.sh "$file"
done