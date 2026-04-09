#!/bin/bash
#SBATCH --output=slurm_logs/%x.%j.out
#SBATCH --error=slurm_logs/%x.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=normal

cd "$SLURM_SUBMIT_DIR"

module purge
module load gaussian/g09-C01

export GAUSS_SCRDIR=/scratch/$USER/$SLURM_JOB_ID
mkdir -p "$GAUSS_SCRDIR"

file="$1"

if [ -z "$file" ]; then
  echo "No input .com file provided"
  exit 1
fi

if [ ! -f "$file" ]; then
  echo "File $file not found"
  exit 1
fi

echo "Running $file"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"

g09 "$file"

echo "End: $(date)"

rm -rf "$GAUSS_SCRDIR"