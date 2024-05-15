#!/bin/bash


#SBATCH --account=SSCM030364
#SBATCH --job-name=simple_job
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-5
#SBATCH --time=03:00:00
#SBATCH --mem=100G

module add languages/r/4.3.1

echo "Run Model Fitting (All except Methylation)"
date
Rscript FeatureSelection_ModelFitting_binary.R
date
echo "End of run"
