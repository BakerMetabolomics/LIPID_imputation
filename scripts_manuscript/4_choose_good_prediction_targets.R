
##############################################################################
## Choose good prediction targets i.e.
## Test predictability of lipids missing in LIPID trial within AusDiab 
## by comparing observed and predicted AusDiab values 
##
## Corresponds to results in Figure 6. and Supplementary Table 3.
##############################################################################

library("tidyverse")
library("magrittr")
library("rsample")
library("glmnet")
library("parallel")
library("doParallel")
library("ggrepel")

##############################################################################
##
##  Load saved working datasets & predictor exclusion lists
##
##############################################################################

## Load saved composite datasets
ausdb_composite <- readRDS("data_derived/Ausdiab_composite_working_data_log_scale.rds")

## Load splits of composite analysis lipid names 
Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')

## Load various predictor exclusion lists from file 2_...; Define matching & missing sets
exclude_predictors_pcor_measure <- readRDS('results/exclude_predictors_pcor_measure.rds')

exclude_pcor_13 <- exclude_predictors_pcor_measure$mad_2$excluded[1:13]

matching_307c_minus <- setdiff(Ausdiab_LIPID_lipid_names_collection$matching_307c, exclude_pcor_13) # 294
missing_522c <- Ausdiab_LIPID_lipid_names_collection$missing_522c



##############################################################################
##
## Within AusDiab assessment of predictions to choose good prediction targets
##
##############################################################################

## Create a 10-fold cross validation split
set.seed(707)
folds <- vfold_cv(data=ausdb_composite, v=10, repeats=1, strata='bmi', breaks=4)


## Function for the job
GlmnetPredictMissingLipids <- function(folds, matching, missing, alpha) {
  
  
  set.seed(707)
  
  ## Loop through folds
  result_ <- lapply(seq(10),function(i_) {
    
    # separate training predictors & targets
    train_x <- analysis(folds$splits[[i_]]) %>% 
      select(id, age, bmi, sex, chol_treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
    train_y <- analysis(folds$splits[[i_]]) %>% 
      select(id, all_of(missing)) %>% column_to_rownames('id') %>% as.matrix()
    
    # separate test predictors & targets 
    test_x <- assessment(folds$splits[[i_]]) %>% 
      select(id, age, bmi, sex, chol_treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
    test_y <- assessment(folds$splits[[i_]]) %>% 
      select(id, all_of(missing)) %>% column_to_rownames('id') %>% as.matrix()
    
    
    ## Loop through lipids, one at a time
    predict_data <- lapply(seq(ncol(train_y)),function(y_) {
      
      ## Build models
      optim_model <- cv.glmnet(x=train_x, y=train_y[,y_,drop=FALSE], family='gaussian', 
                               alpha=alpha, 
                               lambda=exp(seq(-7,1.25,length.out=100)), standardize=TRUE, parallel=TRUE)
      
      
      ## Use optimal lambda to predict lipid concentrations
      y_prediction <- predict(optim_model, test_x, s="lambda.min")
      
      ## Extract betas (nozeros)
      coefs <- as.matrix(coef(optim_model, s="lambda.min"))
      
      ## Get model statistics
      stats <- c(lipid=colnames(train_y)[y_], 
                 fold=i_, 
                 pearson=cor(y_prediction,test_y[,y_]), 
                 coefs[,1])
      
      ## Return predicted data
      return(list('pred'=y_prediction, 'stats'=stats))
      
    })
    
    ## Merge statistics
    stats <- exec(rbind, !!!map(predict_data, 'stats'))
    
    ## Merge predicted lipid data
    predicted_y <- exec(cbind, !!!map(predict_data, 'pred'))
    
    ## Update column names
    colnames(predicted_y) <- colnames(test_y)
    
    ## Return data
    return(list('pred'=predicted_y, 'stats'=stats))
    
  })
  

  ## Get all statistics
  all_stats <- exec(rbind, !!!map(result_, 'stats')) # 5280 x 308
  all_stats <- as_tibble(all_stats) %>% mutate(across(.cols=2:302, .fns=as.numeric))
  
  
  ## Stack all predictions, reorder them as in dat
  preds_y <- exec(rbind, !!!map(result_, 'pred')) # 9922 x 528
  preds_y <- preds_y[as.character(ausdb_composite$id),] 
  
  
  ## Get overall correlations
  corr <- diag(cor(x=preds_y, 
                   y=select(ausdb_composite, id, all_of(missing)) %>% 
                     column_to_rownames('id') %>% as.matrix()))
  
  return(list('corr'=corr, 'stats'=all_stats))
  
}



result_a0.0 <- GlmnetPredictMissingLipids(folds, matching_307c_minus, missing_522c, alpha=0.0)
result_a0.1 <- GlmnetPredictMissingLipids(folds, matching_307c_minus, missing_522c, alpha=0.1)
result_a0.25 <- GlmnetPredictMissingLipids(folds, matching_307c_minus, missing_522c, alpha=0.25)
result_a0.5 <- GlmnetPredictMissingLipids(folds, matching_307c_minus, missing_522c, alpha=0.5)
result_a0.75 <- GlmnetPredictMissingLipids(folds, matching_307c_minus, missing_522c, alpha=0.75)

############################################
format(object.size(result_a0.0), units = "MB")


saveRDS(result_a0.0, 'results/predictability_ausdb_cv_result_a0.0.rds') 
saveRDS(result_a0.1, 'results/predictability_ausdb_cv_result_a0.1.rds') 
saveRDS(result_a0.25, 'results/predictability_ausdb_cv_result_a0.25.rds') 
saveRDS(result_a0.5, 'results/predictability_ausdb_cv_result_a0.5.rds') 
saveRDS(result_a0.75, 'results/predictability_ausdb_cv_result_a0.75.rds') 



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



## Combine all & asses ------------------------------------------------------##
corr522_alpha_series <- tibble(lipid=names(result_a0.0$corr), 
                               a_0=result_a0.0$corr, 
                               a_0.1=result_a0.1$corr, 
                               a_0.25=result_a0.25$corr, 
                               a_0.5=result_a0.5$corr, 
                               a_0.75=result_a0.75$corr)

corr522_alpha_series <- corr522_alpha_series %>% 
  rowwise() %>% 
  mutate(max_corr=names(cur_data()[-1])[which.max(c_across(a_0:a_0.75))]) %>% 
  mutate(best_alpha=as.numeric(str_sub(max_corr, 3,-1)))

corr522_alpha_series <- corr522_alpha_series %>% 
  rowwise() %>% 
  mutate(max=max(c_across(a_0:a_0.75))) %>% 
  rename(best_col=max_corr) %>% 
  rename(max_corr=max) %>% 
  relocate(max_corr, .after='a_0.75')


corr522_alpha_series %>% filter(max_corr>0.7) %>% nrow() # 0.5: 512 (-10; 98%); 0.6: 496 (-26; 95%); 0.7: 475 (-47; 91%)


write_csv(corr522_alpha_series, 'results/predict_ausdb_corr522_alpha_series.csv')
# corr522_alpha_series <- read_csv('results/predict_ausdb_corr522_alpha_series.csv')


## The plan is to predict/impute in LIPID all of those with AusDiab corr > 0.58 (522-20 = 502)
## then it is easy to exclude at more stringent measure (0.6 or 0.7) later
predictable_lipids_alphas <- corr522_alpha_series %>% 
  filter(max_corr > 0.58) %>% 
  select(lipid, max_corr, best_alpha)

count(predictable_lipids_alphas, best_alpha)

# 1       0       59
# 2       0.1    213
# 3       0.25   131
# 4       0.5     58
# 5       0.75    41



## Basic plots ----------------------------------------------------------------
##-----------------------------------------------------------------------------

## Histogram of correlations
ggplot(corr522_alpha_series, aes(max_corr)) +
  geom_histogram(binwidth=0.0125, colour="white", size=0.1, fill="firebrick2", alpha=0.9) + 
  geom_vline(xintercept=0.6, colour="dodgerblue", size=0.25) + 
  theme_bw(base_size=7) + 
  scale_x_continuous(breaks=seq(0,1,0.1)) + 
  labs(x="Correlation")
ggsave("figures/adj_histogram_predict_ausdb_corr522_6.pdf", width=9, height=6, scale=1, units="cm")



# list worst 48 species (corr < 0.7)
corr522_alpha_series %>% filter(max_corr<=0.7) %>% arrange(max_corr) %>% select(lipid, max_corr, best_alpha) %>% print(n=50)

## Table excluded prediction targets
excluded_prediction_targets <- corr522_alpha_series %>% 
  filter(max_corr<=0.6) %>% 
  arrange(max_corr) %>% 
  mutate(max_corr=round(max_corr,3)) %>% 
  select(lipid, max_corr, best_alpha) %>% 
  rename('Lipid'=lipid, 'Max. correlation'=max_corr, 'Best alpha'=best_alpha)

write.xlsx(excluded_prediction_targets, file = "results/excluded_targets_table.xlsx")

##-----------------------------------------------------------------------------

