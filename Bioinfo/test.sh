#!/bin/bash
#SBATCH --job-name=test_openmm_cuda
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=00:10:00
#SBATCH --output=test_openmm_cuda.out
#SBATCH --error=test_openmm_cuda.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate openff

nvidia-smi
python test_openmm_cuda.py