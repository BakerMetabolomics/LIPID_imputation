
##############################################################################
##
##  Multiple Prediction-Imputation of lipids MEASURED in LIPID using Ausdiab
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



##############################################################################
##
##  General function for the prediction/imputation of measured LIPID lipids
##
##############################################################################

GlmnetPredictLipidUsingAusdiab <- function(seed, ausdb.data, lipid.data, matching) {
  
  se7 <- 7*(10*seed)^2 + 7
  set.seed(se7)
  
  # set predictors and targets
  ausdb_x <- ausdb.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  
  lipid_x <- lipid.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  
  # set containers
  lipid_y <- tibble(id=rownames(lipid_x), .rows=nrow(lipid_x))
  lipid_y_sto <- tibble(id=rownames(lipid_x), .rows=nrow(lipid_x))
  ausdb_corr <- vector("numeric", length=ncol(ausdb_x))
  
  # loop through lipids
  for (y_ in seq(ncol(ausdb_x))) {
    
    y_name <- colnames(ausdb_x)[y_]
    
    # build model
    cv.fit <- cv.glmnet(x=ausdb_x[,-y_], 
                        y=ausdb_x[,y_,drop=FALSE], 
                        family='gaussian', 
                        alpha=0.1, 
                        nfolds=50, 
                        lambda=exp(seq(-7,1.25,length.out=250)), 
                        standardize=TRUE, 
                        parallel=TRUE, 
                        trace.it=0)
    
    # predict
    y_pred <- predict(cv.fit, lipid_x[,-y_], s="lambda.1se")
    lipid_y[y_name] <- y_pred
    
    # predict stochastic using internal cv error variance
    rn <- rnorm(n=length(y_pred), mean=0, sd=sqrt(cv.fit$cvm[cv.fit$index[1]]))
    y_pred_rn <- y_pred + matrix(rn, ncol=1)
    
    lipid_y_sto[y_name] <- y_pred_rn
    
    # check accuracy of ausdb predictions just in case
    y_hat <- predict(cv.fit, ausdb_x[,-y_], s="lambda.1se")
    ausdb_corr[y_] <- cor(y_hat, ausdb_x[,y_,drop=FALSE])
    
  }
  
  # remove age, sex, bmi, treat
  lipid_y <- lipid_y %>% select(-c(age, sex, bmi, treat))
  lipid_y_sto <- lipid_y_sto %>% select(-c(age, sex, bmi, treat))
  ausdb_corr <- ausdb_corr[-(1:4)]
  
  return(list('LIPID_miss'=lipid_y, 'LIPID_miss_sto'=lipid_y_sto, 'corr_check'=ausdb_corr))
}



##############################################################################
##############################################################################
##
##  Final measured LIPID multi prediction/imputation
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
                           matching=matching_307c_minus)



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


saveRDS(multi_LIPID_predict, 'results/multi_LIPID_predicted_measured_species_294.rds')
# multi_LIPID_predict <- readRDS('results/multi_LIPID_predicted_measured_species_294.rds')


multi_LIPID_predict[[1]]$LIPID_miss
multi_LIPID_predict[[1]]$LIPID_miss_sto


print(object.size(multi_LIPID_predict), units='auto')

