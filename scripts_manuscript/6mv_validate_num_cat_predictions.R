
##############################################################################
##
##  Validate with multi imputed missing LIPID species: Predict num / cat vars
##
##  Prediction to be done with:
##  1. measured LIPID species
##  2. imputed (unmeasured) LIPID species only
##  3. measured + imputed species
##  4. multiple imputed (unmeasured) LIPID species only
##  5. measured + multiple imputed species
##
##  Corresponds to results in Table 1.
##############################################################################

library("tidyverse")
library("magrittr")
library("rsample")
library("glmnet")
library("parallel")
library("doParallel")
library("pROC")
library("openxlsx")


##############################################################################
##
##  Load saved multiple imputed missing LIPID species
##
##############################################################################

multi_LIPID_predict <- readRDS('results/multi_LIPID_predicted_species_502.rds')

multi_LIPID_predict[[1]]$LIPID_miss
multi_LIPID_predict[[1]]$LIPID_miss_sto


##############################################################################
##
##  Load saved working datasets & predictor, target exclusion lists
##
##############################################################################

## Load saved composite datasets
lipid_composite <- readRDS("data_derived/LIPID_composite_working_data_log_scale.rds")

## Harmonise names of treatment variables
lipid_composite <- lipid_composite %>% rename(treat=treat_12)

## Load splits of composite analysis lipid names 
Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')


## 294 predictors
##---------------------------------------------------------------------------##
## Chose predictors to exclude (from analysis in file 2_..)
exclude_predictors_pcor_measure <- readRDS('results/exclude_predictors_pcor_measure.rds')
exclude_pcor_13 <- exclude_predictors_pcor_measure$mad_2$excluded[1:13]

matching_307c_minus <- setdiff(Ausdiab_LIPID_lipid_names_collection$matching_307c, exclude_pcor_13) # 294


## 414 targets (out of 502)
##---------------------------------------------------------------------------##
## Chose predictable targets (from analysis in file 4_..)
corr522_alpha_series <- read_csv('results/predict_ausdb_corr522_alpha_series.csv')

predictable_lipids_alphas <- corr522_alpha_series %>% 
  filter(max_corr > 0.6) %>% 
  pull(lipid)



##############################################################################
##############################################################################
##
##  Create LIPID data versions (wt/wo imputed) for validation predictions 
##
##############################################################################
##############################################################################



##############################################################################
##  Prepare original LIPID trial data
##############################################################################

## Note: few lipid names bellow will need to be changed to be in line with matching AusDiab names from 12/2021
lipid <- read_csv(file="/Volumes/labs/Metabolomics/Projects/LIPID/Databases/2017_09_12 - Lipid Trial Database - Mapped to Ausdiab.csv", 
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf) %>% 
  rename(`PC(39:5) (b)`=`PC (39:5) (b)`, `PE(O-18:0/22:5)`=`PE(O-18:0/22:5) (a)`, `PC(28:0)`=`PC 28:0`) %>% 
  rename(sbp=SBP, chol=chol0q)

baseline_id <- lipid %>% filter(period==0) %>% mutate(id=as.character(ID)) %>% pull(id)
followup_id <- lipid %>% filter(period==12) %>% mutate(id=as.character(ID)) %>% pull(id)

## select only baseline, exclude 13 discordant lipids and add cvd variable
lipid_composite_base_cvd <- left_join(
  lipid_composite %>% 
    filter(id %in% baseline_id) %>% 
    select(-all_of(exclude_pcor_13)), 
  lipid %>% 
    mutate(id=as.character(ID)) %>% 
    filter(id %in% baseline_id) %>% 
    select(id, sbp, chol, eventCVD, eventStroke), 
  by="id") %>% 
  relocate(sbp, chol, eventCVD, eventStroke, .after=treat)



###############################################################################
## Test: Prepare LIPID + imputed data (exact predictions)
###############################################################################

## Be aware that 82 "non-composite matching" are now doubled in imputations and in original LIPID data
## Remove them from predicted set - They will be taken from LIPID as "closest measured match"

## Overall exclude these from predicted species
lipid_names <- lipid %>% select(`CE(14:0)`:`Hex3Cer(d18:1/24:1)`) %>% colnames() # 342
exclude_from_predicted <- c(intersect(lipid_names, colnames(multi_LIPID_predict[[1]]$LIPID_miss)))

lipid_orig_imputed <- left_join(
  lipid_composite_base_cvd, 
  multi_LIPID_predict[[1]]$LIPID_miss %>% 
    select(id, all_of(predictable_lipids_alphas)) %>% 
    select(-all_of(exclude_from_predicted)) %>% 
    mutate(across(.cols=`AC(12:0)`:`Ubiquinone`, .fns=scale)), 
  by="id"
) # to iterate this w/i a function below


##############################################################################
##  Test: Prepare imputed lipids only data (exact predictions)
##############################################################################

lipid_imputed_only <- lipid_orig_imputed %>% 
  select(-all_of(matching_307c_minus)) # to iterate this w/i a function below


## To avoid complication of setting this up again, make lists of 
## desired original LIPID lipids (294) and imputed LIPID lipids (414):
matching_307c_minus # 294
predictable_414_lipids <- lipid_imputed_only %>% 
  select(-c(id, age, sex, bmi, treat, sbp, chol, eventCVD, eventStroke)) %>% 
  colnames() # 414




##############################################################################
##############################################################################
##
##  Validation for NUMERICAL variables 
##
##############################################################################
##############################################################################


###############################################################################
## Function for numerical variables
###############################################################################
GlmnetPredictCVNum <- function(folds, target, alpha) {
  
  fold_corr <- vector("numeric", length=nrow(folds))
  fold_r2 <- vector("numeric", length=nrow(folds))
  
  all_predictions <- vector("numeric", length=0)
  all_observed <- vector("numeric", length=0)
  
  predictors <- analysis(folds$splits[[1]]) %>% select(-c(id,treat,sbp,chol,eventCVD,eventStroke)) %>% colnames()
  
  # set.seed(707)
  
  ## Loop through folds
  for (i_ in seq(nrow(folds)))  {
    
    # separate training predictors & targets
    train_x <- analysis(folds$splits[[i_]]) %>% 
      select(all_of(predictors)) %>% as.matrix()
    train_y <- analysis(folds$splits[[i_]]) %>% 
      select(all_of(target)) %>% as.matrix()
    
    # separate test predictors & targets 
    test_x <- assessment(folds$splits[[i_]]) %>% 
      select(all_of(predictors)) %>% as.matrix()
    test_y <- assessment(folds$splits[[i_]]) %>% 
      select(all_of(target)) %>% as.matrix()
    
    
    ## Fit the model
    optim_model <- cv.glmnet(x=train_x, y=train_y[,1,drop=FALSE], family='gaussian', 
                             alpha=alpha, 
                             nlambda=200, standardize=TRUE, parallel=TRUE)
    
    ## Predict
    prediction <- predict(optim_model, newx=test_x, s="lambda.min")[,1]
    observed <- test_y[,1]
    
    ## Append and calculate fold correlation 
    all_predictions <- append(all_predictions, prediction)
    all_observed <- append(all_observed, observed)
    
    fold_corr[i_] <- cor(prediction, observed)
    fold_r2[i_] <- (cor(prediction, observed))^2
    
  }
  
  ## Overal 
  oll_corr <- cor(all_predictions, all_observed)
  oll_r2 <- (cor(all_predictions, all_observed))^2
  
  corr_se <- sd(fold_corr, na.rm=T) / sqrt(length(fold_corr))
  r2_se <- sd(fold_r2, na.rm=T) / sqrt(length(fold_r2))
  
  oll_stats <- c(oll_corr, corr_se, oll_r2, r2_se)
  names(oll_stats) <- c("Correlation", "cor_se", "R2", "R2_se")
  
  return(list('overall_stats'=oll_stats, 'fold_corr'=fold_corr, 'fold_r2'=fold_r2, 
              'fold_predictions'=all_predictions, 'fold_observed'=all_observed))
  
}



###############################################################################
## Function to do the job on 10 x 5 data versions 
###############################################################################
GlmnetMultiPredictCVNum <- function(orig_data, multi_imputed_data, predictable_lipids, target_var) {
  
  ## set variable selectors
  lipid_names <- orig_data %>% select(`CE(14:0)`:`TG(58:8) [NL-22:6]`) %>% colnames()
  
  ## set container
  result_5x_10 <- list()
  
  ## loop
  for(j_ in seq(length(multi_imputed_data))) {
    
    # set LIPID + imputed data
    orig_imputed <- left_join(orig_data, 
                              multi_imputed_data[[j_]]$LIPID_miss %>% 
                                select(id, all_of(predictable_lipids)) %>% 
                                mutate(across(.cols=`AC(12:0)`:`Ubiquinone`, .fns=scale)), 
                              by="id")
    
    # set imputed data only
    imputed_only <- orig_imputed %>% 
      select(-all_of(lipid_names))
    
    # set LIPID + stochastic imputed data
    orig_imputed_sto <- left_join(orig_data, 
                                  multi_imputed_data[[j_]]$LIPID_miss_sto %>% 
                                    select(id, all_of(predictable_lipids)) %>% 
                                    mutate(across(.cols=`AC(12:0)`:`Ubiquinone`, .fns=scale)), 
                                  by="id")
    
    # set stochastic imputed data only
    imputed_only_sto <- orig_imputed_sto %>% 
      select(-all_of(lipid_names))
    
    if (j_ == 1) {
      
      # bind all in a list
      data_5x <- list('original'=orig_data,
                      'imputed'=imputed_only, 'original_imputed'=orig_imputed,
                      'imputed_sto'=imputed_only_sto, 'original_imputed_sto'=orig_imputed_sto)
      
      # make folds
      set.seed(707)
      folds_5x <- map(.x=data_5x, .f=vfold_cv, v=10, repeats=1, strata='bmi', breaks=4)
      
      # do the analysis ("original" was unnecessarily repeated every time..)
      result <- map(.x=folds_5x, .f=GlmnetPredictCVNum, target=target_var, alpha=0.1)
      
    } else {
      
      # bind all in a list
      data_2x <- list('imputed_sto'=imputed_only_sto,'original_imputed_sto'=orig_imputed_sto)
      
      # make folds
      set.seed(707)
      folds_2x <- map(.x=data_2x, .f=vfold_cv, v=10, repeats=1, strata='bmi', breaks=4)
      
      # do the analysis ("original" was unnecessarily repeated every time..)
      result <- map(.x=folds_2x, .f=GlmnetPredictCVNum, target=target_var, alpha=0.1)
      
    }
    
    result_5x_10[[j_]] <- result
    
  }
  
  return(result_5x_10)
  
}



###############################################################################
##      Make functions to summary the results above      ##
###############################################################################

# Helper functions to fisher transform r etc. -----------##
fisher_transform <- function(x){
  z <- 0.5 * log((1+x)/(1-x))
  return(z)
}

se_sq <- function(x){
  y <- var(x)/length(x)
  return(y)
}

fisher_back_transform <- function(x){
  r <- (exp(2*x)-1)/(exp(2*x)+1)
  return(r)
}

# Summary function (simplify further... later) ----------##
SummaryMultiPredictCVNum <- function(result){
  
  # "Obtaining Predictions from Models Fit to Multiply Imputed Data"
  # Andrew Miles, Sociological Methods & Research (2016)
  
  # set container
  result_summary <- matrix(0,5,4)
  colnames(result_summary) <- c("r", "lo", "hi", "SE")
  rownames(result_summary) <- c("original","imputed","original+imputed","mult.imputed","original+mult.imputed")
  
  # original -----------------------------------#
  r_j_ave <- result[[1]]$original$overall_stats["Correlation"]
  z_j_ave <- fisher_transform(r_j_ave)
  
  r_ij <- result[[1]]$original$fold_corr
  z_ij <- fisher_transform(r_ij)
  
  se_w <- sd(z_ij)/sqrt(length(z_ij)) # across 10 cv folds
  
  r_dn <- fisher_back_transform(z_j_ave - qnorm(p=0.975)*se_w)
  r_up <- fisher_back_transform(z_j_ave + qnorm(p=0.975)*se_w) 
  
  result_summary[1,] <- c(r_j_ave, r_dn, r_up, se_w)
  
  # imputed ------------------------------------#
  r_j_ave <- result[[1]]$imputed$overall_stats["Correlation"]
  z_j_ave <- fisher_transform(r_j_ave)
  
  r_ij <- result[[1]]$imputed$fold_corr
  z_ij <- fisher_transform(r_ij)
  
  se_w <- sd(z_ij)/sqrt(length(z_ij)) # across 10 cv folds
  
  r_dn <- fisher_back_transform(z_j_ave - qnorm(p=0.975)*se_w)
  r_up <- fisher_back_transform(z_j_ave + qnorm(p=0.975)*se_w) 
  
  result_summary[2,] <- c(r_j_ave, r_dn, r_up, se_w)
  
  # original+imputed ---------------------------#
  r_j_ave <- result[[1]]$original_imputed$overall_stats["Correlation"]
  z_j_ave <- fisher_transform(r_j_ave)
  
  r_ij <- result[[1]]$original_imputed$fold_corr
  z_ij <- fisher_transform(r_ij)
  
  se_w <- sd(z_ij)/sqrt(length(z_ij)) # across 10 cv folds
  
  r_dn <- fisher_back_transform(z_j_ave - qnorm(p=0.975)*se_w)
  r_up <- fisher_back_transform(z_j_ave + qnorm(p=0.975)*se_w) 
  
  result_summary[3,] <- c(r_j_ave, r_dn, r_up, se_w)
  
  # mult.imputed -------------------------------#
  r_j <- map(result, 'imputed_sto') %>% map('overall_stats') %>% map('Correlation') %>% unlist()
  z_j <- fisher_transform(r_j)
  z_j_ave <- mean(z_j)
  
  r_ij <- map(result, 'imputed_sto') %>% map('fold_corr')
  z_ij <- r_ij %>% map(fisher_transform)
  
  # within imputations variance (across 10 cv folds)
  var_w <- z_ij %>% map(se_sq) %>%unlist() %>% mean()
  
  # between imputations variance
  var_b <- sum((z_j-z_j_ave)^2) / (length(z_j)-1)
  
  # total variance and se
  var_oll <- var_w + var_b*(1+(1/length(z_j)))
  se_oll <- sqrt(var_oll)
  
  # r and 95% CI on back-transformed scale
  r_j_ave <- fisher_back_transform(z_j_ave)
  r_dn <- fisher_back_transform(z_j_ave - qnorm(p=0.975)*se_oll)
  r_up <- fisher_back_transform(z_j_ave + qnorm(p=0.975)*se_oll)
  
  result_summary[4,] <- c(r_j_ave, r_dn, r_up, se_oll)
  
  # original+mult.imputed ----------------------#
  r_j <- map(result, 'original_imputed_sto') %>% map('overall_stats') %>% map('Correlation') %>% unlist()
  z_j <- fisher_transform(r_j)
  z_j_ave <- mean(z_j)
  
  r_ij <- map(result, 'original_imputed_sto') %>% map('fold_corr')
  z_ij <- r_ij %>% map(fisher_transform)
  
  # within imputations variance
  var_w <- z_ij %>% map(se_sq) %>%unlist() %>% mean()
  
  # between imputations variance
  var_b <- sum((z_j-z_j_ave)^2) / (length(z_j)-1)
  
  # total variance and se
  var_oll <- var_w + var_b*(1+(1/length(z_j)))
  se_oll <- sqrt(var_oll)
  
  # r and 95% CI on back-transformed scale
  r_j_ave <- fisher_back_transform(z_j_ave)
  r_dn <- fisher_back_transform(z_j_ave - qnorm(p=0.975)*se_oll)
  r_up <- fisher_back_transform(z_j_ave + qnorm(p=0.975)*se_oll)
  
  result_summary[5,] <- c(r_j_ave, r_dn, r_up, se_oll)
  
  ## get summary table
  result_summary <- result_summary %>% as.data.frame() %>% rownames_to_column("data set") %>% as_tibble()
  return(result_summary)
  
}



###############################################################################
## Run the prediction function on few numerical outcomes

## Set clusters -----------------------------------------##
detectCores(all.tests = FALSE, logical = TRUE)
cluster <- makeCluster(6, type='FORK')
registerDoParallel(cl=cluster, cores=6)
getDoParWorkers() 


## SBP
result_sbp <- GlmnetMultiPredictCVNum(orig_data=lipid_composite_base_cvd, 
                                      multi_imputed_data=multi_LIPID_predict, 
                                      predictable_lipids=predictable_414_lipids, 
                                      target_var="sbp")

## Chol
result_chol <- GlmnetMultiPredictCVNum(orig_data=lipid_composite_base_cvd, 
                                       multi_imputed_data=multi_LIPID_predict, 
                                       predictable_lipids=predictable_414_lipids, 
                                       target_var="chol")
# these are saved below along with validations for categorical variables (cvd, stroke)



###############################################################################
## Run the summary function to extract results

## Summary - SBP
summary_sbp <- SummaryMultiPredictCVNum(result_sbp)
summary_sbp %>% mutate('ci.spread'=hi-lo)
sbp <- summary_sbp %>% 
  mutate(across(2:4, .fns=round, digits=3)) %>% 
  mutate(across(5, .fns=round, digits=4)) %>% 
  unite(col="CI", lo:hi, sep="-") %>% 
  rename(Correlation=r) %>% rename("Data set"='data set')

## Summary - Chol
summary_chol <- SummaryMultiPredictCVNum(result_chol)
summary_chol %>% mutate('ci.spread'=hi-lo)
chol <- summary_chol %>% 
  mutate(across(2:4, .fns=round, digits=3)) %>% 
  mutate(across(5, .fns=round, digits=4)) %>% 
  unite(col="CI", lo:hi, sep="-") %>% 
  rename(Correlation=r) %>% rename("Data set"='data set')
# these are saved below along with summaries for categorical variables (cvd, stroke)




##############################################################################
##############################################################################
##
##  Validation for CATEGORICAL variables 
##
##############################################################################
##############################################################################


###############################################################################
## Function for categorical variables
###############################################################################
GlmnetPredictCVCat <- function(folds, target, measure, alpha) {
  
  fold_auc <- vector("numeric", length=nrow(folds))
  
  all_observed <- vector("numeric", length=0)
  all_predprobs <- vector("numeric", length=0)
  all_predictions <- vector("numeric", length=0)
  
  predictors <- analysis(folds$splits[[1]]) %>% select(-c(id,treat,sbp,chol,eventCVD,eventStroke)) %>% colnames()
  
  # set.seed(707)
  
  ## Loop through folds
  for (i_ in seq(nrow(folds)))  {
    
    # separate training predictors & targets
    train_x <- analysis(folds$splits[[i_]]) %>% 
      select(all_of(predictors)) %>% as.matrix()
    train_y <- analysis(folds$splits[[i_]]) %>% 
      select(all_of(target)) %>% as.matrix()
    
    # separate test predictors & targets 
    test_x <- assessment(folds$splits[[i_]]) %>% 
      select(all_of(predictors)) %>% as.matrix()
    test_y <- assessment(folds$splits[[i_]]) %>% 
      select(all_of(target)) %>% as.matrix()
    
    
    ## Fit the model
    optim_model <- cv.glmnet(x=train_x, y=train_y[,1,drop=FALSE], family='binomial', type.measure=measure, 
                             alpha=alpha, 
                             nlambda=200, standardize=TRUE, parallel=TRUE)
    
    ## Predict
    predprob <- predict(optim_model, newx=test_x, s="lambda.min", type="response")[,1]
    observed <- test_y[,1]
    
    ## Get ROC object
    roc_obj <- roc(response=observed, predictor=predprob, quiet=TRUE)
    fold_auc[i_] <- roc_obj$auc[1]
    
    ## Append and calculate fold correlation 
    all_observed <- append(all_observed, observed)
    all_predprobs <- append(all_predprobs, predprob)
    
  }
  
  ## Overall 
  oll_roc_obj <- roc(response=all_observed, predictor=all_predprobs)
  
  oll_coords_y <- coords(oll_roc_obj, x="best", input="threshold", 
                         ret=c("threshold","specificity","sensitivity","npv", "ppv", "tn","tp","fn","fp"), 
                         best.method="youden", best.weights=c(1, 0.5))
  auc_se <- sd(fold_auc, na.rm=T)/sqrt(length(fold_auc))
  oll_stats <- c(matrix(oll_roc_obj$auc[1], dimnames=list(NULL,"auc"))[1,], 
                 matrix(auc_se, dimnames=list(NULL,"auc_se"))[1,], 
                 as.matrix(oll_coords_y)[1,])
  
  oll_prediction <- if_else(all_predprobs >= oll_stats["threshold"], 1, 0)
  
  
  return(list('overall_stats'=oll_stats, 'roc'=oll_roc_obj, 'fold_auc'=fold_auc, 
              'predictions'=oll_prediction, 'observed'=all_observed))
  
}



###############################################################################
## Function to do the job on 10 x 5 data versions 
###############################################################################
GlmnetMultiPredictCVCat <- function(orig_data, multi_imputed_data, predictable_lipids, target_var) {
  
  ## set variable selectors
  lipid_names <- orig_data %>% select(`CE(14:0)`:`TG(58:8) [NL-22:6]`) %>% colnames()
  
  ## set container
  result_5x_10 <- list()
  
  ## loop
  for(j_ in seq(length(multi_imputed_data))) {
    
    # set LIPID + imputed data
    orig_imputed <- left_join(orig_data, 
                              multi_imputed_data[[j_]]$LIPID_miss %>% 
                                select(id, all_of(predictable_lipids)) %>% 
                                mutate(across(.cols=`AC(12:0)`:`Ubiquinone`, .fns=scale)), 
                              by="id")
    
    # set imputed data only
    imputed_only <- orig_imputed %>% 
      select(-all_of(lipid_names))
    
    # set LIPID + stochastic imputed data
    orig_imputed_sto <- left_join(orig_data, 
                                  multi_imputed_data[[j_]]$LIPID_miss_sto %>% 
                                    select(id, all_of(predictable_lipids)) %>% 
                                    mutate(across(.cols=`AC(12:0)`:`Ubiquinone`, .fns=scale)), 
                                  by="id")
    
    # set stochastic imputed data only
    imputed_only_sto <- orig_imputed_sto %>% 
      select(-all_of(lipid_names))
    
    if (j_ == 1) {
      
      # bind all in a list
      data_5x <- list('original'=orig_data,
                      'imputed'=imputed_only, 'original_imputed'=orig_imputed,
                      'imputed_sto'=imputed_only_sto, 'original_imputed_sto'=orig_imputed_sto)
      
      # make folds
      set.seed(707)
      folds_5x <- map(.x=data_5x, .f=vfold_cv, v=10, repeats=1, strata='bmi', breaks=4)
      
      # do the analysis ("original" was unnecessarily repeated every time..)
      result <- map(.x=folds_5x, .f=GlmnetPredictCVCat, target=target_var, measure='auc', alpha=0.1)
      
    } else {
      
      # bind all in a list
      data_2x <- list('imputed_sto'=imputed_only_sto,'original_imputed_sto'=orig_imputed_sto)
      
      # make folds
      set.seed(707)
      folds_2x <- map(.x=data_2x, .f=vfold_cv, v=10, repeats=1, strata='bmi', breaks=4)
      
      # do the analysis ("original" was unnecessarily repeated every time..)
      result <- map(.x=folds_2x, .f=GlmnetPredictCVCat, target=target_var, measure='auc', alpha=0.1)
      
    }
    
    result_5x_10[[j_]] <- result
    
  }
  
  return(result_5x_10)
  
}



###############################################################################
##      Make functions to summary the results above      ##
###############################################################################

# Helper functions etc. -----------##
se_sq <- function(x){
  y <- var(x)/length(x)
  return(y)
}



# Summary function (simplify further... later) ----------##
SummaryMultiPredictCVCat <- function(result){
  
  # "Obtaining Predictions from Models Fit to Multiply Imputed Data"
  # Andrew Miles, Sociological Methods & Research (2016)
  
  # set container
  result_summary <- matrix(0,5,4)
  colnames(result_summary) <- c("AUC", "lo", "hi", "SE")
  rownames(result_summary) <- c("original","imputed","original+imputed","mult.imputed","original+mult.imputed")
  
  # original -----------------------------------#
  auc_j_ave <- result[[1]]$original$overall_stats["auc"]
  
  se_w <- result[[1]]$original$overall_stats["auc_se"] # across 10 cv folds
  
  r_dn <- auc_j_ave - qnorm(p=0.975)*se_w
  r_up <- auc_j_ave + qnorm(p=0.975)*se_w
  
  result_summary[1,] <- c(auc_j_ave, r_dn, r_up, se_w)
  
  # imputed ------------------------------------#
  auc_j_ave <- result[[1]]$imputed$overall_stats["auc"]
  
  se_w <- result[[1]]$imputed$overall_stats["auc_se"] # across 10 cv folds
  
  r_dn <- auc_j_ave - qnorm(p=0.975)*se_w
  r_up <- auc_j_ave + qnorm(p=0.975)*se_w
  
  result_summary[2,] <- c(auc_j_ave, r_dn, r_up, se_w)
  
  # original+imputed ---------------------------#
  auc_j_ave <- result[[1]]$original_imputed$overall_stats["auc"]
  
  se_w <- result[[1]]$original_imputed$overall_stats["auc_se"] # across 10 cv folds
  
  r_dn <- auc_j_ave - qnorm(p=0.975)*se_w
  r_up <- auc_j_ave + qnorm(p=0.975)*se_w
  
  result_summary[3,] <- c(auc_j_ave, r_dn, r_up, se_w)
  
  # mult.imputed -------------------------------#
  auc_j <- map(result, 'imputed_sto') %>% map('overall_stats') %>% map('auc') %>% unlist()
  auc_j_ave <- mean(auc_j)
  
  auc_ij <- map(result, 'imputed_sto') %>% map('fold_auc')
  
  # within imputations variance (across 10 cv folds)
  var_w <- auc_ij %>% map(se_sq) %>%unlist() %>% mean()
  
  # between imputations variance
  var_b <- sum((auc_j-auc_j_ave)^2) / (length(auc_j)-1)
  
  # total variance and se
  var_oll <- var_w + var_b*(1+(1/length(auc_j)))
  se_oll <- sqrt(var_oll)
  
  # r and 95% CI on back-transformed scale
  r_dn <- auc_j_ave - qnorm(p=0.975)*se_oll
  r_up <- auc_j_ave + qnorm(p=0.975)*se_oll
  
  result_summary[4,] <- c(auc_j_ave, r_dn, r_up, se_oll)
  
  # original+mult.imputed ----------------------#
  auc_j <- map(result, 'original_imputed_sto') %>% map('overall_stats') %>% map('auc') %>% unlist()
  auc_j_ave <- mean(auc_j)
  
  auc_ij <- map(result, 'original_imputed_sto') %>% map('fold_auc')
  
  # within imputations variance (across 10 cv folds)
  var_w <- auc_ij %>% map(se_sq) %>%unlist() %>% mean()
  
  # between imputations variance
  var_b <- sum((auc_j-auc_j_ave)^2) / (length(auc_j)-1)
  
  # total variance and se
  var_oll <- var_w + var_b*(1+(1/length(auc_j)))
  se_oll <- sqrt(var_oll)
  
  # r and 95% CI on back-transformed scale
  r_dn <- auc_j_ave - qnorm(p=0.975)*se_oll
  r_up <- auc_j_ave + qnorm(p=0.975)*se_oll
  
  result_summary[5,] <- c(auc_j_ave, r_dn, r_up, se_oll)
  
  
  ## get summary table
  result_summary <- result_summary %>% as.data.frame() %>% rownames_to_column("data set") %>% as_tibble()
  return(result_summary)
  
}



###############################################################################
## Run the prediction function on few categorical outcomes

## Set clusters -----------------------------------------##
detectCores(all.tests = FALSE, logical = TRUE)
cluster <- makeCluster(6, type='FORK')
registerDoParallel(cl=cluster, cores=6)
getDoParWorkers() 


## CVD
result_cvd <- GlmnetMultiPredictCVCat(orig_data=lipid_composite_base_cvd, 
                                      multi_imputed_data=multi_LIPID_predict, 
                                      predictable_lipids=predictable_414_lipids, 
                                      target_var="eventCVD")

## Stroke
result_stroke <- GlmnetMultiPredictCVCat(orig_data=lipid_composite_base_cvd, 
                                         multi_imputed_data=multi_LIPID_predict, 
                                         predictable_lipids=predictable_414_lipids, 
                                         target_var="eventStroke")


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


## Save all validation results together -------------------------------------##
validation_multi_impute_results <- list('SBP'=result_sbp, 'Cholesterol'=result_chol, 
                                        'CVD'=result_cvd, 'Stroke'=result_stroke)
saveRDS(validation_multi_impute_results, 'results/validation_multi_impute_results_variation.rds')
# validation_multi_impute_results <- readRDS('results/validation_multi_impute_results_variation.rds')



###############################################################################
## Run the summary function to extract results

## Summary - CVD
summary_cvd <- SummaryMultiPredictCVCat(result_cvd)
summary_cvd %>% mutate('ci.spread'=hi-lo)
cvd <- summary_cvd %>% 
  mutate(across(2:4, .fns=round, digits=3)) %>% 
  mutate(across(5, .fns=round, digits=4)) %>% 
  unite(col="CI", lo:hi, sep="-") %>% 
  rename("Data set"='data set')

## Summary - SBtroke
summary_stroke <- SummaryMultiPredictCVCat(result_stroke)
summary_stroke %>% mutate('ci.spread'=hi-lo)
stroke <- summary_stroke %>% 
  mutate(across(2:4, .fns=round, digits=3)) %>% 
  mutate(across(5, .fns=round, digits=4)) %>% 
  unite(col="CI", lo:hi, sep="-") %>% 
  rename("Data set"='data set')




###############################################################################
##  Write the result for all 4 variables
###############################################################################

## List of validation summary tables ----------------------------------------##
validation_summ_tables <- list('SBP'=sbp, 
                               'Cholesterol'=chol, 
                               'CVD'=cvd, 
                               'Stroke'=stroke)

## Write a multi-sheet excel file
write.xlsx(validation_summ_tables, file="results/validation_multi_impute_summ_tables_variation.xlsx")


