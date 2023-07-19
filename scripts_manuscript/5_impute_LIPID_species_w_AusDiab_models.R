
##############################################################################
##
##  Prediction-Imputation of lipids missing in LIPID using Composite Ausdiab
##
##############################################################################

library("tidyverse")
library("magrittr")
library("glmnet")
library("parallel")
library("doParallel")

##############################################################################
##
##  Load saved working datasets & predictor, target exclusion lists
##
##############################################################################

## Load saved composite datasets
ausdb_composite <- readRDS("data_derived/Ausdiab_composite_working_data_log_scale.rds")
lipid_composite <- readRDS("data_derived/LIPID_composite_working_data_log_scale.rds")

## Harmonise names of treatment variables
ausdb_composite <- ausdb_composite %>% rename(treat=chol_treat)
lipid_composite <- lipid_composite %>% rename(treat=treat_12)


## Load splits of composite analysis lipid names 
Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')


## 294 predictors
##---------------------------------------------------------------------------##
## Chose predictors to exclude (from analysis in file 2_..)
exclude_predictors_pcor_measure <- readRDS('results/exclude_predictors_pcor_measure.rds')
exclude_pcor_13 <- exclude_predictors_pcor_measure$mad_2$excluded[1:13]

matching_307c_minus <- setdiff(Ausdiab_LIPID_lipid_names_collection$matching_307c, exclude_pcor_13) # 294


## 496 targets
##---------------------------------------------------------------------------##
## Chose predictable targets and their optimal alpha (from analysis in file 4_..)
corr522_alpha_series <- read_csv('results/predict_ausdb_corr522_alpha_series.csv')

predictable_lipids_alphas <- corr522_alpha_series %>% 
  filter(max_corr > 0.6) %>% 
  select(lipid, max_corr, best_alpha)

missing_alpha_522c_minus <- predictable_lipids_alphas$best_alpha
names(missing_alpha_522c_minus) <- predictable_lipids_alphas$lipid

## The aim is to predict/impute in LIPID all of those with AusDiab corr > 0.6



##############################################################################
##
##  General function for the final prediction/imputation
##
##############################################################################

GlmnetPredictLipidUsingAusdiab <- function(ausdb.data, lipid.data, matching, missing, alpha.vec) {
  
  # ausdb.data = reference data used to build the prediction models 
  # lipid.data = target data to predict lipids into 
  # matching = character vector of predictor names (lipids) that were sufficiently consistent between reference and target data
  # missing = character vector of target names (lipids) to be predicted that were well-predicted when tested in reference data
  # alpha.vec = numerical vector of pre-tested best performing elastic net alpha values for each individual target lipid (alt. set all to 0.1 or so)
  
  # set predictors and targets
  ausdb_x <- ausdb.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  ausdb_y <- ausdb.data %>% select(id, all_of(missing)) %>% column_to_rownames('id') %>% as.matrix()
  
  lipid_x <- lipid.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  
  # set containers
  lipid_y <- tibble(id=rownames(lipid_x), .rows=nrow(lipid_x))
  lipid_y_ave <- tibble(id=rownames(lipid_x), .rows=nrow(lipid_x))
  ausdb_corr <- vector("numeric", length=ncol(ausdb_y))
  
  # loop through lipids
  for (y_ in seq(ncol(ausdb_y))) {
    
    y_name <- colnames(ausdb_y)[y_]
    alpha <- unname(alpha.vec[y_name])
    
    # build model
    cv.fit <- cv.glmnet(x=ausdb_x, 
                        y=ausdb_y[,y_,drop=FALSE], 
                        family='gaussian', 
                        alpha=alpha, 
                        nfolds=10, 
                        lambda=exp(seq(-7,1.25,length.out=300)), 
                        standardize=TRUE, 
                        parallel=TRUE, 
                        trace.it=0)
    
    # predict
    y_pred <- predict(cv.fit, lipid_x, s="lambda.1se")
    lipid_y[y_name] <- y_pred
    
    # predict naive stochastic x5 using internal cv error variance
    rn <- rnorm(n=length(y_pred)*5, mean=0, sd=sqrt(cv.fit$cvm[cv.fit$index[1]]))
    y_pred_rn <- matrix(rep(y_pred, times=5), ncol=5) + matrix(rn, ncol=5)
    
    y_pred_ave <- matrix(rowMeans(y_pred_rn), ncol=1)
    lipid_y_ave[y_name] <- y_pred_ave
    
    # check accuracy of ausdb predictions just in case
    y_hat <- predict(cv.fit, ausdb_x, s="lambda.1se")
    ausdb_corr[y_] <- cor(y_hat, ausdb_y[,y_,drop=FALSE])
    
  }
  
  return(list('LIPID_miss'=lipid_y, 'LIPID_miss_ave'=lipid_y_ave, 'corr_check'=ausdb_corr))
}



##############################################################################
##############################################################################
##
##  Final LIPID prediction/imputation
##
##############################################################################
##############################################################################

## Set clusters -----------------------------------------##
detectCores(all.tests = FALSE, logical = TRUE)
cluster <- makeCluster(6, type='FORK')
registerDoParallel(cl=cluster, cores=6)
getDoParWorkers() 


## Run predictions 
set.seed(707)

LIPID_predict <- GlmnetPredictLipidUsingAusdiab(ausdb.data=ausdb_composite, 
                                                lipid.data=lipid_composite, 
                                                matching=matching_307c_minus, 
                                                missing=names(missing_alpha_522c_minus), 
                                                alpha.vec=missing_alpha_522c_minus)


## Stop clusters ----------------------------------------##
stopCluster(cl=cluster)
stopImplicitCluster()
getDoParWorkers() 

# in case of parallel trouble run this:
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()


saveRDS(LIPID_predict, 'results/LIPID_predicted_species_496.rds')
# LIPID_predict <- readRDS('results/LIPID_predicted_species_496.rds')




##############################################################################
##############################################################################
##
##  Create a new LIPID data with predicted/imputed lipid species
##
##############################################################################
##############################################################################


##############################################################################
##  Load original LIPID trial & AusDiab
##############################################################################

## Note: few lipid names bellow will need to be changed to be in line with matching AusDiab names from 12/2021
lipid <- read_csv(file="/Volumes/labs/Metabolomics/Projects/LIPID/Databases/2017_09_12 - Lipid Trial Database - Mapped to Ausdiab.csv", 
                       na = c("", "NA"), trim_ws = TRUE, guess_max = Inf)

## Fix some lipid names and select only matching species (for the purpose of prediction)
lipid <- lipid %>% 
  rename(`PC(39:5) (b)`=`PC (39:5) (b)`, `PE(O-18:0/22:5)`=`PE(O-18:0/22:5) (a)`, `PC(28:0)`=`PC 28:0`)

## Some checks
sum(Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82 %in% colnames(lipid)) ## 82
all.equal(intersect(Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82, colnames(lipid)), 
          intersect(Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82, Ausdiab_LIPID_lipid_names_collection$matching_307i)) # TRUE

## AusDiab 
ausdb <- read_csv(file="/Volumes/labs/Metabolomics/Projects/Ausdiab/Databases/2021_12_16 - Ausdiab lipidomics and clinical data combined.csv", 
                  col_select = c(id, drage_00, drsex_00_n, bmi_00, choltabl_00, q20_angi_00, q20_coro_00, q20_stro_00, CVDdnfev_10, cvd_d, cva_d, 
                                 `AC(12:0)`:`Ubiquinone`), 
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf)


##############################################################################
##  Prepare variables
##############################################################################


## Be aware that 82 "non-composite matching" are now doubled in imputations and in original LIPID data
## Remove them from predicted set - They will be taken from LIPID as "closest measured match"
## BTW, apart from PC(O-40:7) (b) they are perfectly correlated i.e. similar
## PC(O-40:7) (b) and PC(O-40:7) (a) need to be removed as their composite version PC(O-40:7) is excluded predictor
intersect(colnames(lipid), colnames(LIPID_predict$LIPID_miss)) # 82
intersect(intersect(colnames(lipid), colnames(LIPID_predict$LIPID_miss)), 
          Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82) # 82


## Also 12 excluded "bad predictors" should be removed from both LIPID and AusDiab when they are used together


## Overall exclude these from original and predicted species ----------------##
exclude_from_predicted <- c(intersect(colnames(lipid), colnames(LIPID_predict$LIPID_miss)), "PC(O-40:7) (a)") # same as: noncomposite_matching_82 + 1 = 83
exclude_from_predictors <- c(setdiff(exclude_pcor_13, c("PC(O-40:7)")), "PC(O-40:7) (b)") # 13; TG(54:0) [NL-18:0] ... do not exclude these, just mark them as AusDb-LIPID discordant



###############################################################################
## Create complete LIPID data with exact predictions of missing species ##
###############################################################################

original_LIPID_names_corrected <- lipid %>% select(`CE(14:0)`:`Hex3Cer(d18:1/24:1)`) %>% colnames() %>% sort()

LIPID_original_lipids <-  
  lipid %>% 
  select(ID, all_of(original_LIPID_names_corrected))
  
LIPID_original_and_imputed_lipids_logged_scaled <- left_join(
  LIPID_original_lipids %>% 
    mutate(across(.cols=`CE(14:0)`:`TG(58:8) [NL-22:6]`, .fns=log)) %>% 
    mutate(across(.cols=`CE(14:0)`:`TG(58:8) [NL-22:6]`, .fns=scale)), 
  LIPID_predict$LIPID_miss %>% 
    select(-all_of(exclude_from_predicted)) %>% 
    mutate(across(.cols=`AC(12:0)`:`Ubiquinone`, .fns=scale)) %>% 
    mutate(id=as.numeric(id)), 
  by=c("ID"="id")
)


## Have a reference of which were original and predicted species ------------##
original_LIPID_species <- original_LIPID_names_corrected # 342

imputed_LIPID_species <- LIPID_original_and_imputed_lipids_logged_scaled %>% 
  select(`AC(12:0)`:`Ubiquinone`) %>% colnames() # 413

expanded_LIPID_species <- union(original_LIPID_species, imputed_LIPID_species) %>% sort() # 755


## Make final for saving; Reorder lipid names
LIPID_original_lipids <- LIPID_original_lipids
  as.data.frame()

LIPID_original_and_imputed_lipids_logged_scaled <- LIPID_original_and_imputed_lipids_logged_scaled %>% 
  select(ID, all_of(expanded_LIPID_species)) %>% 
  as.data.frame()


## Save data sets and reference of lipid type (original vs. imputed etc.) ---##
write.csv(LIPID_original_lipids, 'data_derived/2023_05_22_LIPID_trial_original_lipids.csv', row.names=FALSE) # 11,773 x 343 (ID+342 lipids)
write.csv(LIPID_original_and_imputed_lipids_logged_scaled, 'data_derived/2023_05_22_LIPID_trial_original_and_imputed_lipids_logged_scaled.csv', row.names=FALSE) # 11,773 x 756 (ID+755 lipids)


## Save reference for original/imputed & AusDiab-LIPID concordant/non-concordant
LIPID_lipid_dictionary <- tibble(lipid=expanded_LIPID_species) %>% 
  mutate(origin=if_else(lipid %in% original_LIPID_species, "original", "imputed")) %>% 
  mutate(AusDiab_discordant=if_else(lipid %in% exclude_from_predictors, "yes", "no")) %>% 
  mutate(AusDiab_absent=if_else(lipid %in% Ausdiab_LIPID_lipid_names_collection$nomap_lipid_names_35, "yes", "no")) %>% 
  as.data.frame()

write.csv(LIPID_lipid_dictionary, 'data_derived/2023_05_22_LIPID_lipid_dictionary.csv', row.names=FALSE) # 755 x 4



###############################################################################
## Add reference of lipid type (original vs. imputed) to lipid names collection
###############################################################################
Ausdiab_LIPID_lipid_names_collection <- 
  c(Ausdiab_LIPID_lipid_names_collection, 
    list('original_LIPID_species' = original_LIPID_species, # 342
         'imputed_LIPID_species' = imputed_LIPID_species, # 413
         'excluded_AusDb_predictors' = c(setdiff(exclude_pcor_13, "PC(O-40:7)"), "PC(O-40:7) (b)"), # 13
         'excluded_AusDb_targets' = corr522_alpha_series %>% filter(max_corr <= 0.6) %>% pull(lipid), # 26
         'congruent_AusDb_LIPID_species' = setdiff(expanded_LIPID_species, 
                                                   union(Ausdiab_LIPID_lipid_names_collection$nomap_lipid_names_35, 
                                                         c(setdiff(exclude_pcor_13, "PC(O-40:7)"), "PC(O-40:7) (b)"))))) # 707


saveRDS(Ausdiab_LIPID_lipid_names_collection, 'data_derived/Ausdiab_LIPID_lipid_names_collection.rds')
# Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')

