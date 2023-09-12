# LIPID Imputation Repository
This repository contains the R scripts used in the manuscript titled "Imputation of plasma lipid species to facilitate integration of lipidomic datasets" by Aleksandar Dakic, Jingqin Wu, Tingting Wang, Kevin Huynh, Natalie Mellett, Thy Duong, Habtamu B Beyene, Dianna J Magliano, Jonathan E Shaw, Melinda J Carrington, Mike Inouye, Jean Y Yang, Gemma A Figtree, John Simes, the LIPID Study Investigators, Corey Giles, and Peter Meikle.

The scripts in this repository demonstrate a framework for harmonising lipidomic datasets with different levels of granularity in their lipid measurements, through imputation.

## Scripts Overview
The scripts are located in the scripts_manuscript folder and are written in R.
- **1_make_working_datasets.R**: Creates the working datasets for the study.
- **2_find_discrepant_species_AusDiab_LIPID.R**: Finds discrepant species in the AusDiab LIPID dataset.
- **3_choose_good_predictors.R**: Selects good predictors for the imputation models.
- **4_choose_good_prediction_targets.R**: Chooses good prediction targets for the imputation models.
- **5_impute_LIPID_species_w_AusDiab_models.R**: Imputes LIPID species with the AusDiab models.
- **5m_multi_impute_measured_LIPID_species.R**: Performs multiple imputations on measured LIPID species.
- **5m_multi_impute_missing_LIPID_species.R**:Performs multiple imputations on missing LIPID species.
- **6m_validate_CVD_lipid_assoc_all.R**: Validates the associations between cardiovascular disease (CVD) and all lipid species.
- **6m_validate_CVD_lipid_assoc_measured.R**: Validates the associations between CVD and measured lipid species.
- **6mv_validate_num_cat_predictions.R**: Validates the numerical and categorical predictions.
- **code_source_data_to figures.R**: Generates figures from the source data.

Simulated data can be produced by running the script **0_create_simulated_datasets.R**, which is inside the simulated_data folder.

## Usage
To use these scripts, clone the repository and run the scripts in the order listed above. Simulated data can be created using the R script in the simulated_data folder.
Please ensure that you have the necessary R packages installed.

## Dependencies
The scripts were run with:
- R (version: 4.1.2)
- R Packages (with versions):
  - broom (v1.0.5)
  - corpcor (v1.6.10)
  - doParallel (v1.0.17)
  - ggrepel (v0.9.3)
  - glmnet (v4.1.7)
  - igraph (v1.4.2)
  - magrittr (v2.0.3)
  - openxlsx (v4.2.5.2)
  - parallel (v4.1.2)
  - pdist (v1.2.1)
  - pROC (v1.18.0)
  - RColorBrewer (v1.1.3)
  - readxl (v1.4.2)
  - rsample (v1.2.0)
  - tableone (v0.13.2)
  - tidyverse (v2.0.0)

These packages are available on CRAN and can be installed using this command:
```
install.packages(c("broom","corpcor","doParallel","ggrepel","glmnet","igraph","magrittr","openxlsx","parallel","pdist","pROC","RColorBrewer","readxl","rsample","tableone","tidyverse"))
```

## Data
The data used in these scripts is not included in this repository. Please refer to the manuscript for information on how to access the individual level data.
