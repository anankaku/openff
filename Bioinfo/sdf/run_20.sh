#!/bin/bash
#SBATCH --job-name=openmm_md
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=openmm_20.out
#SBATCH --error=openmm_20.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate openff

python run_20.py