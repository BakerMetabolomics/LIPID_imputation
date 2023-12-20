
## Sex-stratified version of Fig.3 analysis (Reviewer 1 suggestion) --> Supp. Figure 2
## (compare with the version where sex is only a predictor)

library("tidyverse")
library("magrittr")
library("glmnet")
library("parallel")
library("doParallel")
library("rsample")


## Load saved composite datasets
ausdb_composite <- readRDS("data_derived/Ausdiab_composite_working_data_log_scale.rds")
lipid_composite <- readRDS("data_derived/LIPID_composite_working_data_log_scale.rds")

## Load splits of composite analysis lipid names  
Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')

## Select variables needed here
ausdb307 <- ausdb_composite %>% 
  select(age, sex, bmi, chol_treat, all_of(Ausdiab_LIPID_lipid_names_collection$matching_307c)) %>% 
  rename(treat=chol_treat)
ausdb307_f <- ausdb307 %>% filter(sex==0)
ausdb307_m <- ausdb307 %>% filter(sex==1)

lipid307 <- lipid_composite %>% 
  select(age, sex, bmi, treat_12, all_of(Ausdiab_LIPID_lipid_names_collection$matching_307c)) %>% 
  rename(treat=treat_12)
lipid307_f <- lipid307 %>% filter(sex==0)
lipid307_m <- lipid307 %>% filter(sex==1)


## Load various predictor exclusion lists from file 2_...
exclude_predictors_pcor_measure <- readRDS('results/exclude_predictors_pcor_measure.rds')

exclude_pcor_13 <- exclude_predictors_pcor_measure$mad_2$excluded[1:13]



### Function for this job ###
#############################

GlmnetPredictIndividualLipids <- function(train.data, test.data, exclude, sex_stratified, alpha) {
  
  predictors <- setdiff(colnames(train.data), exclude)
  
  ### Mixed sex models, but segregate prediction accuracy
  if(sex_stratified=="no") {
    train <- train.data %>% select(all_of(predictors)) %>% as.matrix()
    test_f <- test.data %>% filter(sex==0) %>% select(all_of(predictors)) %>% as.matrix()
    test_m <- test.data %>% filter(sex==1) %>% select(all_of(predictors)) %>% as.matrix()
    
    ## Loop through lipids, one at a time - Females
    predict_data_f <- sapply(seq(ncol(train)),function(y_) {
      
      ## Build ridge model
      optim_model <- cv.glmnet(x=train[,-y_], 
                               y=train[,y_,drop=FALSE], 
                               family='gaussian', 
                               alpha=alpha, 
                               lambda=exp(seq(-7,1.25,length.out=200)), 
                               standardize=TRUE, 
                               parallel=TRUE)
      
      ## Use optimal lambda to predict lipid concentrations
        y_prediction_f <- predict(optim_model, test_f[,-y_], s="lambda.1se")
   
      ## Get model statistics
      pearson <- cor(y_prediction_f, test_f[,y_])
      return(pearson)
      
    })
    names(predict_data_f) <- predictors
    
    
    ## Loop through lipids, one at a time - Males
    predict_data_m <- sapply(seq(ncol(train)),function(y_) {
      
      ## Build ridge model
      optim_model <- cv.glmnet(x=train[,-y_], 
                               y=train[,y_,drop=FALSE], 
                               family='gaussian', 
                               alpha=alpha, 
                               lambda=exp(seq(-7,1.25,length.out=200)), 
                               standardize=TRUE, 
                               parallel=TRUE)
      
      ## Use optimal lambda to predict lipid concentrations
      y_prediction_m <- predict(optim_model, test_m[,-y_], s="lambda.1se")
      
      ## Get model statistics
      pearson <- cor(y_prediction_m, test_m[,y_])
      return(pearson)
      
    })
    names(predict_data_m) <- predictors
    
    return(list(female=predict_data_f, male=predict_data_m))
    
  
  ### Sex-stratified models, and prediction accuracy  
  } else {
    train <- train.data %>% select(all_of(predictors)) %>% select(-sex) %>% as.matrix()
    test <- test.data %>% select(all_of(predictors)) %>% select(-sex) %>% as.matrix()
    
    ## Loop through lipids, one at a time
    predict_data <- sapply(seq(ncol(train)),function(y_) {
      
      ## Build ridge model
      optim_model <- cv.glmnet(x=train[,-y_], 
                               y=train[,y_,drop=FALSE], 
                               family='gaussian', 
                               alpha=alpha, 
                               lambda=exp(seq(-7,1.25,length.out=200)), 
                               standardize=TRUE, 
                               parallel=TRUE)
      
      ## Use optimal lambda to predict lipid concentrations
      y_prediction <- predict(optim_model, test[,-y_], s="lambda.1se")
      
      ## Get model statistics
      pearson <- cor(y_prediction, test[,y_])
      return(pearson)
      
    })
    names(predict_data) <- setdiff(predictors, "sex")
    
    return(predict_data)
  }
  
} # just ignore first 3-4 predictions as they are for age, sex, bmi, treatment



## Set clusters ----------------------------------------##
detectCores(all.tests = FALSE, logical = TRUE)
cluster <- makeCluster(6, type='FORK')
registerDoParallel(cl=cluster, cores=6)
getDoParWorkers()



##############################################################################
## Only run the analysis on LIPID predictions (mixed sex vs. sex-stratified) 
## Only run the analysis with 13 bad predictors removed
##############################################################################

## Run LIPID mixed sex analysis ---------------------------------------------##
set.seed(707)
res_lipid <- GlmnetPredictIndividualLipids(train.data=ausdb307, 
                                                   test.data=lipid307, 
                                                   exclude=exclude_pcor_13, 
                                                   sex_stratified="no", 
                                                   alpha=0.1)
# extract predicted-observed correlations
hist(res_lipid$female[-c(1,2,3,4)], breaks=100)
hist(res_lipid$male[-c(1,2,3,4)], breaks=100)


## Run LIPID female-specific analysis ---------------------------------------##
set.seed(707)
res_lipid_f <- GlmnetPredictIndividualLipids(train.data=ausdb307_f, 
                                                   test.data=lipid307_f, 
                                                   exclude=exclude_pcor_13, 
                                                   sex_stratified="yes", 
                                                   alpha=0.1)
# extract predicted-observed correlations
hist(res_lipid_f[-c(1,2,3)], breaks=100)


## Run LIPID male-specific analysis -----------------------------------------##
set.seed(707)
res_lipid_m <- GlmnetPredictIndividualLipids(train.data=ausdb307_m, 
                                             test.data=lipid307_m, 
                                             exclude=exclude_pcor_13, 
                                             sex_stratified="yes", 
                                             alpha=0.1)
# extract predicted-observed correlations
hist(res_lipid_m[-c(1,2,3)], breaks=100)



## Stop clusters  ---------------------------------------##
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()



## GGPlot & save the results ------------------------------------------------##
sex_data <- tibble(mixed_f=res_lipid$female[-c(1,2,3,4)],
                   strat_f=res_lipid_f[-c(1,2,3)], 
                   mixed_m=res_lipid$male[-c(1,2,3,4)], 
                   strat_m=res_lipid_m[-c(1,2,3)])

# female
ggplot(data=sex_data, 
       aes(x=mixed_f, y=strat_f)) + 
  geom_point(size=1.5, alpha=0.4, colour="dodgerblue", stroke=0.4) + 
  geom_abline(slope=1, intercept=0, colour="grey30", linewidth=0.3) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(limits=c(0.25,1), breaks=seq(0.2,1,0.2)) + 
  scale_y_continuous(limits=c(0.25,1), breaks=seq(0.2,1,0.2)) + 
  labs(x="mixed-sex model", y="sex-stratified model", title="Female") + 
  theme(axis.text=element_text(size=5), plot.title=element_text(hjust=0.5))
ggsave("figures/sex_stratified_F_supp2.pdf", width=6.25, height=6.5, scale=1, units="cm")

# male
ggplot(data=sex_data, 
       aes(x=mixed_m, y=strat_m)) + 
  geom_point(size=1.5, alpha=0.4, colour="tomato", stroke=0.4) + 
  geom_abline(slope=1, intercept=0, colour="grey30", linewidth=0.3) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(limits=c(0.25,1), breaks=seq(0.2,1,0.2)) + 
  scale_y_continuous(limits=c(0.25,1), breaks=seq(0.2,1,0.2)) + 
  labs(x="mixed-sex model", y="sex-stratified model", title="Male") + 
  theme(axis.text=element_text(size=5), plot.title=element_text(hjust=0.5))
ggsave("figures/sex_stratified_M_supp2.pdf", width=6.25, height=6.5, scale=1, units="cm")



