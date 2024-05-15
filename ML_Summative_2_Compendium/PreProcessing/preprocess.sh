#!/bin/bash


#SBATCH --account=SSCM030364
#SBATCH --job-name=simple_job
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=50G

module add languages/r/4.3.1

echo "Run Preprocessing"
date
Rscript PreProcessing.R
date
echo "End of run"


