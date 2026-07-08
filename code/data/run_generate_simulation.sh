#!/bin/bash
#SBATCH --job-name=sim_app_generate
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sim_app_generate_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
##SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name


set -euo pipefail

# Make sure the log directory exists (SLURM needs it to exist before the job starts writing)
mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs

# Activate the conda environment used by sim_app
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sim-app-tutorial

cd /dcs04/hicks/data/Jan/sim_project/sim_paper

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: $(which python)"

python generate_simulation.py

echo "Job finished: $(date)"
