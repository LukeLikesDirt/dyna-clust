#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# Exicute the R script
Rscript ./02.Check_infraspecific_annotations.R

# Deactivate conda environment