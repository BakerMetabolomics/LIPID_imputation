
##############################################################################
## Perform imputations and assesments of accuracy (AusDiab --> SanAntonio)
## 
## Test transferability of AusDiab models to SAFHS by predicting the same set 
## of masked lipids in AusDiab and SAFHS, using the models built in AusDiab 
##
## Corresponds to results in Figure 9. and Supp. Figure 4.
##############################################################################

library("BakerMetabolomics")
library("tidyverse")
library("magrittr")
library("pdist")
library("corpcor")
library("glmnet")
library("parallel")
library("doParallel")
library("rsample")
library("ggrepel")



##############################################################################
##
##  Load saved working datasets & predictor exclusion lists
##
##############################################################################

## Load & prepare AusDb data
loadAusdiab()
Ausdiab <- Ausdiab %>% rownames_to_column("ID") %>% as_tibble() %>% select(-starts_with("Total"))
Ausdiab_covar <- Ausdiab_covar %>% as_tibble()

ausdb <- Ausdiab_covar %>% 
  select(id, .AGE, .SEX, .BMI, choltabl_00) %>% 
  mutate(id=as.character(id)) %>% 
  rename("ID"="id") %>% 
  filter(!choltabl_00=="Missing") %>% 
  mutate(lipid_treat=if_else(choltabl_00=="Yes", 1, 0)) %>% 
  select(ID, .AGE, .SEX, .BMI, lipid_treat) %>% 
  left_join(Ausdiab, by="ID") %>% 
  mutate(across(.cols=!c(ID, .AGE, .SEX, .BMI, lipid_treat), .fns=log)) %>% 
  mutate(across(.cols=!c(ID, .AGE, .SEX, .BMI, lipid_treat), .fns=scale)) # can add more covariates later if change plan 

## Load & prepare SAFHS data
loadSAFHS()
SAFHS <- SAFHS %>% rownames_to_column("ID_Date") %>% as_tibble() %>% select(-starts_with("Total"))
SAFHS_covar <- SAFHS_covar %>% as_tibble()

safhs <- SAFHS_covar %>% 
  select(ID_Date, ID, Visit, Cohort, .AGE, .SEX) %>% 
  left_join(SAFHS, by="ID_Date") %>% 
  mutate(across(.cols=!c(ID_Date, ID, Visit, Cohort, .AGE, .SEX), .fns=log)) %>% 
  mutate(across(.cols=!c(ID_Date, ID, Visit, Cohort, .AGE, .SEX), .fns=scale)) %>% 
  filter(!is.na(.AGE)) %>% filter(!is.na(.SEX)) # also, consider to filter age to min 25 if required


## Note: analysis done on lipid species without missing values and observations without missing data for sex and age
## BMI and lipid treatment were not used as predictors as they were missing on large number of observations

safhs_missing_cols <- sapply(safhs, function(x) which(is.na(x))) # 62 lipids missing 2090-3500 obs

# make overall complete data i.e. remove the above lipids with many missing observations
safhs_complete <- safhs %>% select(which(colSums(is.na(.))==0)) # n=5590, 795 lipids



##############################################################################
##
##  Chose which lipids to blind and impute: 1) same/similar set as for LIPID 
##
##############################################################################

## How many overlapping lipids in total
common_lipids <- intersect(colnames(ausdb)[6:751], colnames(safhs_complete)[7:801]) # 700

## Load lipids names collection from Ausdiab-LIPID work
Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')

# some checks
intersect(Ausdiab_LIPID_lipid_names_collection$imputed_LIPID_species, 
          common_lipids) %>% length() # 377 of 413 with good prediction accuracy in AusDiab (cor>=0.6)
setdiff(Ausdiab_LIPID_lipid_names_collection$imputed_LIPID_species, 
        common_lipids) # 36 mismatch, but some of them due to name updates; manually fix them: 

## Set targets for prediction & predictor lipids (names) --------------------##
target_lipids <- c(intersect(Ausdiab_LIPID_lipid_names_collection$imputed_LIPID_species, 
                             common_lipids), 
                   "LPC(15-MHDA) [sn1] & LPC(17:0) [sn2]", 
                   "LPC(18:3) [sn1] (a) & LPC(18:3) [sn2] (b)", 
                   "LPC(19:0) [sn1] (a) & LPC(19:0) [sn2] (b)", 
                   "PI(15-MHDA_18:1) & PI(17:0_18:1)", 
                   "PI(15-MHDA_18:2) & PI(17:0_18:2)", 
                   "PI(15-MHDA_20:4) & PI(17:0_20:4)", 
                   "SM(d18:1/22:0) & SM(d16:1/24:0)") %>% sort() # 384

predictor_lipids <- setdiff(common_lipids, target_lipids) %>% sort() # 316



##############################################################################
##
## Find irreconcilable species between AusDb & SAFHS to remove from predictors 
##
##############################################################################

###############################################################################
##
## Remove discrepant predictor lipids between two datasets one-by-one based on 
## differences in their Partial correlation vectors. Run pcor after every removal, 
## remove the most discrepant lipid, stop when there's no apparent outliers
##
###############################################################################

## Use a function for this task from script 2_find_discrepant_species...-----##

# input for the function are matrices of same lipids

## Run the function 
lipids_to_exclude_1.75mad <- RemoveDissimilarLipidsPcor(ausdb %>% select(all_of(predictor_lipids)) %>% as.matrix(), 
                                                        safhs_complete %>% select(all_of(predictor_lipids)) %>% as.matrix(), 
                                                        mad, 1.75) # 25
lipids_to_exclude_2sd <- RemoveDissimilarLipidsPcor(ausdb %>% select(all_of(predictor_lipids)) %>% as.matrix(), 
                                                    safhs_complete %>% select(all_of(predictor_lipids)) %>% as.matrix(), 
                                                    sd, 2) # 19

exclude_predictors_pcor_measure <- list('mad_1.75'=lipids_to_exclude_1.75mad, 'sd_2'=lipids_to_exclude_2sd)


#-----------------------------------------------------------------------------#
## Diagnostic plots

## Plot overall Frobenius distances -----------------------------------------##
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
ggplot(tibble(x=0:20, y=exclude_predictors_pcor_measure$mad_1.75$frobenius_ave[1:21]), aes(x, y)) +
  geom_line(size=0.4) + 
  geom_point(size=1.75, shape=21, fill="firebrick1", colour="white", stroke=0.75) +
  theme_bw(base_size=8) + 
  scale_y_continuous(label=scientific_10) + 
  scale_x_continuous(breaks=seq(0,25,5)) + 
  labs(x="Number of excluded lipids", y="Average squared distance") + 
  geom_vline(aes(xintercept=6), linetype="dotted", size=0.4) + 
  geom_vline(aes(xintercept=12), linetype="dotted", size=0.4) + 
  theme(axis.text=element_text(size=6))   
ggsave("figures/adj_average_pcor_sq_dist_ausdb_safhs_supp4b.pdf", width=8, height=6, scale=1, units="cm")

scales::scientific
scale_y_continuous(label=function(x) format(x, scientific = TRUE))

## Plot distribution of per lipid distances for the first 15 exclusions (lipids_to_exclude_1.75mad)
dist_df <- tibble(dist=unlist(exclude_predictors_pcor_measure$mad_1.75$pcor_dist_series[1:16])^2, 
                  excluded=rep(0:15, times=316:(316-15))) %>% 
  mutate(excluded=factor(as.character(excluded), levels=as.character(0:15)))

ggplot(dist_df, aes(dist)) +
  geom_histogram(binwidth=0.03, colour="white", size=0.003, fill="royalblue4") + 
  theme_bw(base_size=7.5) + 
  labs(x="Squared distance") + 
  scale_y_continuous(breaks=c(0,25,50)) + 
  facet_wrap(~excluded, scales="fixed", ncol=4) + 
  theme(axis.text=element_text(size=5)) 
ggsave("figures/adj_distribution_of_pcor_sq_dist_ausdb_safhs_supp4a.pdf", width=10, height=8, scale=1, units="cm")


## 12 predictors to exclude
excluded_predictors <- exclude_predictors_pcor_measure$mad_1.75$excluded[1:12]

## Update predictor lipids !
##############################################################################
predictor_lipids <- setdiff(predictor_lipids, excluded_predictors) %>% sort() # 304 !



##############################################################################
##  Function for the job 
##############################################################################
GlmnetPredictAndAssessTargetLipids <- function(reference, target, predictors, targets, alpha) {
  
  # split train, test for reference out of sample test
  reference_split <- initial_split(reference, prop=0.8, strata=".AGE", breaks=4)
  reference_train <- training(reference_split) 
  reference_test  <- testing(reference_split) 
  
  # set data blocks:
  ## whole reference
  reference_predictors <- reference %>% select(ID, .AGE, .SEX, all_of(predictors)) %>% 
    column_to_rownames('ID') %>% as.matrix()
  reference_targets <- reference %>% select(ID, all_of(targets)) %>% 
    column_to_rownames('ID') %>% as.matrix()
  
  ## whole target
  target_predictors <- target %>% select(ID_Date, .AGE, .SEX, all_of(predictors)) %>% 
    column_to_rownames('ID_Date') %>% as.matrix()
  target_targets <- target %>% select(ID_Date, all_of(targets)) %>% 
    column_to_rownames('ID_Date') %>% as.matrix()
  
  ## reference train split
  reference_predictors_train <- reference_train %>% select(ID, .AGE, .SEX, all_of(predictors)) %>% 
    column_to_rownames('ID') %>% as.matrix()
  reference_targets_train <- reference_train %>% select(ID, all_of(targets)) %>% 
    column_to_rownames('ID') %>% as.matrix()
  ## reference test split
  reference_predictors_test <- reference_test %>% select(ID, .AGE, .SEX, all_of(predictors)) %>% 
    column_to_rownames('ID') %>% as.matrix()
  reference_targets_test <- reference_test %>% select(ID, all_of(targets)) %>% 
    column_to_rownames('ID') %>% as.matrix()
  
  # set containers
  predicted_reference_targets <- tibble(ID=rownames(reference_targets), .rows=nrow(reference_targets))
  predicted_target_targets <- tibble(ID_Date=rownames(target_targets), .rows=nrow(target_targets))
  predicted_reference_targets_test <- tibble(ID=rownames(reference_targets_test), .rows=nrow(reference_targets_test))
  
  # loop through target lipids
  for (y_ in seq(ncol(reference_targets))) {
    
    y_name <- colnames(reference_targets)[y_]
    
    # build model
    cv.fit <- cv.glmnet(x=reference_predictors, 
                        y=reference_targets[,y_,drop=FALSE], 
                        family='gaussian', 
                        alpha=alpha, 
                        nfolds=10, 
                        nlambda=300, 
                        standardize=TRUE, 
                        parallel=TRUE, 
                        trace.it=0)
    
    # predict for reference data
    y_pred_reference <- predict(cv.fit, reference_predictors, s="lambda.1se")
    predicted_reference_targets[y_name] <- y_pred_reference
    
    # predict for target data
    y_pred_target <- predict(cv.fit, target_predictors, s="lambda.1se")
    predicted_target_targets[y_name] <- y_pred_target
    
  }
  
  # loop through target lipids: out of sample reference prediction version
  for (y_ in seq(ncol(reference_targets_train))) {
    
    y_name <- colnames(reference_targets_train)[y_]
    
    # build model
    cv.fit <- cv.glmnet(x=reference_predictors_train, 
                        y=reference_targets_train[,y_,drop=FALSE], 
                        family='gaussian', 
                        alpha=alpha, 
                        nfolds=10, 
                        nlambda=300, 
                        standardize=TRUE, 
                        parallel=TRUE, 
                        trace.it=0)
    
    # predict for reference data
    y_pred_reference_test <- predict(cv.fit, reference_predictors_test, s="lambda.1se")
    predicted_reference_targets_test[y_name] <- y_pred_reference_test
    
  }
  
  # get correlations between observed & predicted
  pred_acc <- tibble(lipid=cor(as_tibble(reference_targets), select(predicted_reference_targets, -ID)) %>% diag() %>% names(), 
                     AusDiab_fit=cor(as_tibble(reference_targets), select(predicted_reference_targets, -ID)) %>% diag(), 
                     AusDiab_test=cor(as_tibble(reference_targets_test), select(predicted_reference_targets_test, -ID)) %>% diag(), 
                     SAFHS=cor(as_tibble(target_targets), select(predicted_target_targets, -ID_Date)) %>% diag())
  
  return(list(prediction_accuracies=pred_acc, 
              AusDb_fit_predictions=predicted_reference_targets, 
              AusDb_test_predictions=predicted_reference_targets_test, 
              SAFHS_predictions=predicted_target_targets))
  
}



##############################################################################
##
##  Predict AusDiab & SAFHS target set; Compare accuracy of predictions 
##
##############################################################################


## Set clusters -----------------------------------------##
detectCores(all.tests = FALSE, logical = TRUE)
cluster <- makeCluster(6, type='FORK')
registerDoParallel(cl=cluster, cores=6)
getDoParWorkers() 


## Run predictions 
set.seed(707)

AusDb_SAFHS_predict <- GlmnetPredictAndAssessTargetLipids(reference=ausdb, 
                                                          target=safhs_complete, 
                                                          predictors=predictor_lipids, 
                                                          targets=target_lipids, 
                                                          alpha=0.1)


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


saveRDS(AusDb_SAFHS_predict, 'results/AusDb_SAFHS_target_predictions_accuracy.rds')



##############################################################################
##
##  Plot the results 
##
##############################################################################


## Figures ----------------------------------------------##
# focused cor 0.4-1
ggplot(data=AusDb_SAFHS_predict$prediction_accuracies, 
       aes(x=AusDiab_test, y=SAFHS)) + 
  geom_point(size=1.75, alpha=0.6, colour="indianred2", stroke=0.4) + 
  geom_label_repel(aes(label=ifelse(abs(AusDiab_test-SAFHS)>0.2, lipid,'')), 
                   box.padding=0.2, label.padding=0.075, point.padding=0.25, label.r=0.1, label.size=0.1, max.time=5, force=2, force_pull=1, 
                   segment.color='grey50', segment.size=0.2, min.segment.length=0.2, size=1.75, nudge_x=-0.01, nudge_y=-0.01, max.overlaps=30) + 
  geom_abline(slope=1, intercept=0, colour="grey30", linewidth=0.4) + 
  geom_abline(slope=1, intercept=-0.2, colour="grey50", linewidth=0.15) + 
  geom_abline(slope=1, intercept=0.2, colour="grey50", linewidth=0.15) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(limits=c(0.4,1), breaks=seq(0.4,1,0.1)) + 
  scale_y_continuous(limits=c(0.4,1), breaks=seq(0.4,1,0.1)) + 
  labs(x="AusDiab", y="SAFHS") + 
  theme(axis.text=element_text(size=6))
ggsave("figures/validation_heterogeneity_Ausdiab12_SAFHS_9a_focused.pdf", width=9, height=9, scale=1, units="cm")

# full view cor 0.1-01
ggplot(data=AusDb_SAFHS_predict$prediction_accuracies, 
       aes(x=AusDiab_test, y=SAFHS)) + 
  geom_point(size=1.75, alpha=0.6, colour="indianred2", stroke=0.4) + 
  geom_label_repel(aes(label=ifelse(abs(AusDiab_test-SAFHS)>0.2, lipid,'')), 
                   box.padding=0.2, label.padding=0.075, point.padding=0.25, label.r=0.1, label.size=0.1, max.time=5, force=2, force_pull=1, 
                   segment.color='grey50', segment.size=0.2, min.segment.length=0.2, size=1.75, nudge_x=-0.01, nudge_y=-0.01, max.overlaps=30) + 
  geom_abline(slope=1, intercept=0, colour="grey30", linewidth=0.4) + 
  geom_abline(slope=1, intercept=-0.2, colour="grey50", linewidth=0.15) + 
  geom_abline(slope=1, intercept=0.2, colour="grey50", linewidth=0.15) + 
  theme_bw(base_size=8) + 
  scale_x_continuous(limits=c(0.08,1), breaks=seq(0.1,1,0.1)) + 
  scale_y_continuous(limits=c(0.08,1), breaks=seq(0.1,1,0.1)) + 
  labs(x="AusDiab", y="SAFHS") + 
  theme(axis.text=element_text(size=6))
ggsave("figures/validation_heterogeneity_Ausdiab12_SAFHS_9a.pdf", width=9, height=9, scale=1, units="cm")

