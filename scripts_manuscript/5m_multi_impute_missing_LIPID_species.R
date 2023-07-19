
##############################################################################
##
##  Multiple Prediction-Imputation of lipids MISSING in LIPID using Ausdiab
##
##############################################################################

library("tidyverse")
library("magrittr")
library("rsample")
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


## 502 targets
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

GlmnetPredictLipidUsingAusdiab <- function(seed, ausdb.data, lipid.data, matching, missing, alpha.vec) {
  
  se7 <- 7*(10*seed)^2 + 7
  set.seed(se7)
  
  # set predictors and targets
  ausdb_x <- ausdb.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  ausdb_y <- ausdb.data %>% select(id, all_of(missing)) %>% column_to_rownames('id') %>% as.matrix()
  
  lipid_x <- lipid.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  
  # set containers
  lipid_y <- tibble(id=rownames(lipid_x), .rows=nrow(lipid_x))
  lipid_y_sto <- tibble(id=rownames(lipid_x), .rows=nrow(lipid_x))
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
                        lambda=exp(seq(-7,1.25,length.out=250)), 
                        standardize=TRUE, 
                        parallel=TRUE, 
                        trace.it=0)
    
    # predict
    y_pred <- predict(cv.fit, lipid_x, s="lambda.1se")
    lipid_y[y_name] <- y_pred
    
    # predict naive stochastic x5 using internal cv error variance
    rn <- rnorm(n=length(y_pred), mean=0, sd=sqrt(cv.fit$cvm[cv.fit$index[1]]))
    y_pred_rn <- y_pred + matrix(rn, ncol=1)
    
    lipid_y_sto[y_name] <- y_pred_rn
    
    # check accuracy of ausdb predictions just in case
    y_hat <- predict(cv.fit, ausdb_x, s="lambda.1se")
    ausdb_corr[y_] <- cor(y_hat, ausdb_y[,y_,drop=FALSE])
    
  }
  
  return(list('LIPID_miss'=lipid_y, 'LIPID_miss_sto'=lipid_y_sto, 'corr_check'=ausdb_corr))
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
multi_LIPID_predict <- map(.x=1:10, .f=GlmnetPredictLipidUsingAusdiab, 
                           ausdb.data=ausdb_composite, 
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


saveRDS(multi_LIPID_predict, 'results/multi_LIPID_predicted_species_502.rds')
# multi_LIPID_predict <- readRDS('results/multi_LIPID_predicted_species_502.rds') 


## If it is easier to work with nested data frame / tibble
multi_LIPID_predict_tib <- left_join(
  enframe(map(multi_LIPID_predict, 'LIPID_miss'), name='index', value='preds_exact'), 
  enframe(map(multi_LIPID_predict, 'LIPID_miss_sto'), name='index', value='preds_stoch')
)


