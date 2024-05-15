#!/bin/bash


#SBATCH --account=SSCM030364
#SBATCH --job-name=simple_job
#SBATCH --partition=teach_cpu
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem=100G

module add languages/r/4.3.1

echo "Run Methylation Model Fitting"
date
Rscript methfinal_analysis.R
date
echo "End of run"


