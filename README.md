
# ML Summative: OMIC and Clinical Data for predicting Cancer Outcomes Repository

**Note:** This repository was restarted multiple times due to the git file corrupting because of large files preventing commits which led to it corrupting. So, there is no log since the start of the project.

## Repository Contents

### `PackagesInstall.R`
- Script for installing required R packages from CRAN; try from GitHub if CRAN installation fails.

### `.gitignore`
- Includes old attempts, mostly very similar codes and plots but before tuning or adjustments to preprocessing.

### `ML_Summative_2_Compendium`:
- Directory containing all code for running primary and exploratory analysis. 

#### `renv' and `renv.lock'
- Contains a record of the computational environment for the project.

#### `rawdata`
- Empty folder for adding raw data, which will then be processed by our R scripts.

#### `preprocessed`
- Repository for saving the preprocessed data. Includes a subfolder for methylation preprocessed sets.

#### `plot`
- Repository to store plots after running all analysis. These include all plots from previous runs.

#### `PreProcessing`
- Contains an R script for preprocessing data from the `rawdata` folder. Accompanied by a .sh file to run this on BC4.
- `methprocess` sub-folder containing R and .sh scripts for methylation preprocessing, which reduce the feature set to 52 features used for prediction.

#### `ModelFitting_withFeatureSelection`
- Contains 2 R scripts and equivalent .sh scripts for running analysis on preprocessed datasets. Separate scripts are used for methylation and other datasets.

#### `PCA`
- Contains an R script and a .sh script to run PCA on combined datasets before model fitting.

## Directions for Running Analysis

1. **Install Required Packages**  
   Run the provided RScript on BC4 to install necessary packages. Or use RENV to replicate computational environemnt.

2. **Copy Raw Data**
   ```bash
   cp path_to_raw_data/*.txt path_to_my_directory/ML_Summative_2_Compendium/rawdata/
   ```

3. **Ensure Required Files Are in `rawdata`**
   Make sure that `protein-annotation.txt`, `protein.txt`, `clinical.txt`, `mrna.txt`, `mutation.txt`, and `methylation.txt` files are all present.

4. **Run Preprocessing Scripts**
   These are in the PreProcessing folder, with some inside methprocess.
   ```bash
   sbatch preprocess.sh
   sbatch meth_processing_job.sh
   sbatch meth_fs_0_job.sh
   # Follow by sbatch for meth_fs_1_job.sh to meth_fs_3_job.sh
   ```

5. **Confirm Preprocessed Data**
   Ensure that the `.Rda` files are present in the `preprocessed` folder.

6. **Run Model Fitting**
   These are in the ModelFitting_withFeatureSelection_andPerformance and running the sh scripts run the corresponding R scripts.
   ```bash
   sbatch FitModels.sh
   sbatch methanalysis_job.sh
   ```

7. **Run PCA**
   Run the PCA script in the PCA folder
   ```bash
   sbatch PCA.sh
   ```

8. **Review Outputs**
   Check the `plot` directory for plots and review slurm files for confusion matrices, accuracies, and F1 scores.


