
##############################################################################
## 
## Test transferability of AusDiab models to LIPID across different 
## lipid-treatment (pravastatin) groups by checking and comparing 
## predictions of lipid species shared by AusDiab & LIPID
##
## Corresponds to results in Figure 4. and Supp. Figure 3.
##############################################################################

library("tidyverse")
library("magrittr")
library("glmnet")
library("parallel")
library("doParallel")
library("rsample")


## Load saved composite datasets
ausdb_composite <- readRDS("data_derived/Ausdiab_composite_working_data_log_scale.rds")
lipid_composite <- readRDS("data_derived/LIPID_composite_working_data_log_scale.rds")

## Harmonise names of treatment variables
ausdb_composite <- ausdb_composite %>% rename(treat=chol_treat)
lipid_composite <- lipid_composite %>% rename(treat=treat_12)

## Load splits of composite analysis lipid names
Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')

## Load original LIPID to extract pid & period (in addition to ID)
lipid <- read_csv(file="data_derived/2017_09_12 - Lipid Trial Database - Mapped to Ausdiab.csv", 
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf) %>% 
  rename(`PC(39:5) (b)`=`PC (39:5) (b)`, `PE(O-18:0/22:5)`=`PE(O-18:0/22:5) (a)`, `PC(28:0)`=`PC 28:0`)


## 294 predictors
##---------------------------------------------------------------------------##
## Chose predictors to exclude (from analysis in script 2_..)
exclude_predictors_pcor_measure <- readRDS('results/exclude_predictors_pcor_measure.rds')
exclude_pcor_13 <- exclude_predictors_pcor_measure$mad_2$excluded[1:13]

matching_307c_minus <- setdiff(Ausdiab_LIPID_lipid_names_collection$matching_307c, exclude_pcor_13) # 294


## Run two versions of analysis using lipids and 4 covariates as predictors (Figure 4)
## vs. only lipids as predictors, i.e. no age, sex, bmi, treatment (Supp. Figure 2):


##############################################################################
##  General function for the prediction/imputation of measured LIPID lipids
##############################################################################

GlmnetPredictLipidUsingAusdiab <- function(seed, ausdb.data, lipid.data, matching) {
  
  se7 <- 7*(10*seed)^2 + 7
  set.seed(se7)
  
  # set predictors and targets
  ausdb_x <- ausdb.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  
  lipid_x <- lipid.data %>% select(id, age, sex, bmi, treat, all_of(matching)) %>% column_to_rownames('id') %>% as.matrix()
  
  # set containers
  lipid_y <- tibble(id=rownames(lipid_x), .rows=nrow(lipid_x))
  ausdb_corr <- vector("numeric", length=ncol(ausdb_x))
  lipid_corr <- vector("numeric", length=ncol(lipid_x))
  
  # loop through lipids
  for (y_ in seq(ncol(ausdb_x))) {
    
    y_name <- colnames(ausdb_x)[y_]
    
    # build model
    cv.fit <- cv.glmnet(x=ausdb_x[,-y_], 
                        y=ausdb_x[,y_,drop=FALSE], 
                        family='gaussian', 
                        alpha=0.1, 
                        nfolds=20, 
                        lambda=exp(seq(-7,1.25,length.out=250)), 
                        standardize=TRUE, 
                        parallel=TRUE, 
                        trace.it=0)
    
    # predict
    y_pred <- predict(cv.fit, lipid_x[,-y_], s="lambda.1se")
    lipid_y[y_name] <- y_pred
    
    
    # check accuracy of ausdb & lipid predictions
    y_hat <- predict(cv.fit, ausdb_x[,-y_], s="lambda.1se")
    ausdb_corr[y_] <- cor(y_hat, ausdb_x[,y_,drop=FALSE])
    lipid_corr[y_] <- cor(y_pred, lipid_x[,y_,drop=FALSE])
    
  }
  
  # remove age, sex, bmi, treat
  lipid_y <- lipid_y %>% select(-c(age, sex, bmi, treat))
  ausdb_corr <- ausdb_corr[-(1:4)]
  lipid_corr <- lipid_corr[-(1:4)]
  
  return(list('LIPID_miss'=lipid_y, 'corr_check_ausdb'=ausdb_corr, 'corr_check_lipid'=lipid_corr))
}



##############################################################################
##############################################################################
##
##  Measured LIPID multi prediction/imputation
##
##############################################################################
##############################################################################

## Set clusters -----------------------------------------##
detectCores(all.tests = FALSE, logical = TRUE)
cluster <- makeCluster(6, type='FORK')
registerDoParallel(cl=cluster, cores=6)
getDoParWorkers() 



## Run predictions 
multi_LIPID_predict <- map(.x=1:2, .f=GlmnetPredictLipidUsingAusdiab, 
                           ausdb.data=ausdb_composite, 
                           lipid.data=lipid_composite, 
                           matching=matching_307c_minus) # with all covars, 2 reps

multi_LIPID_predict_nocovars <- map(.x=1:2, .f=GlmnetPredictLipidUsingAusdiab, 
                                    ausdb.data=ausdb_composite, 
                                    lipid.data=lipid_composite, 
                                    matching=matching_307c_minus) # no covars, 2 reps (adjust the function above)



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


## Save
saveRDS(multi_LIPID_predict, 'results/LIPID_predicted_measured_species_294.rds')
saveRDS(multi_LIPID_predict_nocovars, 'results/LIPID_predicted_nocovars_measured_species_294.rds')

multi_LIPID_predict <- readRDS('results/LIPID_predicted_measured_species_294.rds')
multi_LIPID_predict <- readRDS('results/LIPID_predicted_nocovars_measured_species_294.rds')

multi_LIPID_predict[[1]]$LIPID_miss


## Include "period" and "pid" variables in results to enable filtering (baseline & follow up)
lipid_composite_paired <- 
  left_join(
    lipid %>% 
      select(ID,filter_paired, period, pid, TREAT) %>% 
      filter(filter_paired==0) %>% 
      mutate(ID=as.character(ID)), 
    lipid_composite, 
    by=join_by(ID==id)
  ) %>% 
  select(ID, pid, period, TREAT, all_of(matching_307c_minus)) %>% 
  arrange(pid, period)

lipid_predicted_paired <- 
  left_join(
    lipid_composite_paired %>% 
      select(ID,pid,period,TREAT), 
    multi_LIPID_predict[[1]]$LIPID_miss %>% 
      select(id, all_of(matching_307c_minus)), 
    by=join_by(ID==id)
  ) %>% 
  arrange(pid, period)


## Get correlations of measured-predicted values (accuracy)--@ baseline and @ followup--stratified on statin treatment
cordata <- tibble(
  baseline_nontreat=map2(lipid_composite_paired %>% filter(period==0 & TREAT==0) %>% select(all_of(matching_307c_minus)), 
                         lipid_predicted_paired %>% filter(period==0 & TREAT==0) %>% select(all_of(matching_307c_minus)), 
                         .f=cor) %>% unlist(), 
  followup_nontreat=map2(lipid_composite_paired %>% filter(period==12 & TREAT==0) %>% select(all_of(matching_307c_minus)), 
                         lipid_predicted_paired %>% filter(period==12 & TREAT==0) %>% select(all_of(matching_307c_minus)), 
                         .f=cor) %>% unlist(), 
  baseline_treat=map2(lipid_composite_paired %>% filter(period==0 & TREAT==1) %>% select(all_of(matching_307c_minus)), 
                      lipid_predicted_paired %>% filter(period==0 & TREAT==1) %>% select(all_of(matching_307c_minus)), 
                      .f=cor) %>% unlist(), 
  followup_treat=map2(lipid_composite_paired %>% filter(period==12 & TREAT==1) %>% select(all_of(matching_307c_minus)), 
                      lipid_predicted_paired %>% filter(period==12 & TREAT==1) %>% select(all_of(matching_307c_minus)), 
                      .f=cor) %>% unlist(), 
  lipid=matching_307c_minus) # names(baseline_all)

write_csv(cordata, 'results/data_fig_4new.csv')


###############################################################################
## Reorganise for plotting placebo vs. treatment across two time points
###############################################################################
cordata_long <- cordata %>% 
  pivot_longer(cols=baseline_nontreat:followup_treat, names_to=c("timepoint", "randomisation"), names_sep="_", values_to="value") %>% 
  pivot_wider(names_from=randomisation, values_from=value) %>% 
  mutate(time=if_else(timepoint=="baseline", "Baseline", "Follow up (only pravastatin arm was exposed to treatment)"))

## Plot
my.colors <- c("royalblue4", "seagreen4")
subset_label <- c(baseline="Baseline", followup="Follow up")
ggplot(data=cordata_long, 
       aes(x=nontreat, y=treat)) + 
  geom_point(size=1.4, alpha=0.5, aes(colour=time), stroke=0.3) + 
  geom_abline(slope=1, intercept=0, colour="grey30") + 
  theme_bw(base_size=7) + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2)) + 
  scale_color_manual(values=my.colors) + 
  labs(x="Placebo randomised", y="Pravastatin randomised") + 
  facet_wrap(~timepoint, scales="fixed", ncol=2, labeller=labeller(timepoint=subset_label)) + 
  theme(strip.background=element_rect(fill="white", colour="white"), strip.text=element_text(face=NULL, size=7), axis.text=element_text(size=5), 
        legend.position="bottom", legend.title=element_blank())
ggsave("figures/statin_comp_acc_timepoints_4a.pdf", width=9, height=6.25, scale=1, units="cm")

# for nocovars version:
my.colors <- c("royalblue1", "seagreen3")
# ...
ggsave("figures/statin_comp_acc_nocovars_timepoints_supp3a.pdf", width=12, height=7, scale=1, units="cm")


###############################################################################
## Reorganise for plotting baseline vs. follow-up across two treatment groups
###############################################################################
cordata_long <- cordata %>% 
  pivot_longer(cols=baseline_nontreat:followup_treat, names_to=c("timepoint", "randomisation"), names_sep="_", values_to="value") %>% 
  pivot_wider(names_from=timepoint, values_from=value) %>% 
  mutate(random=if_else(randomisation=="nontreat", "Placebo", "Pravastatin (only follow up time point exposed to treatment)"))

## Plot
my.colors <- c("royalblue4", "seagreen4")
subset_label <- c(nontreat="Placebo randomised", treat="Pravastatin randomised")
ggplot(data=cordata_long, 
       aes(x=baseline, y=followup)) + 
  geom_point(size=1.4, alpha=0.5, aes(colour=random), stroke=0.3) + 
  geom_abline(slope=1, intercept=0, colour="grey30") + 
  theme_bw(base_size=7) + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2)) + 
  scale_color_manual(values=my.colors) + 
  labs(x="Baseline", y="Follow up") + 
  facet_wrap(~randomisation, scales="fixed", ncol=3, labeller=labeller(randomisation=subset_label)) + 
  theme(strip.background=element_rect(fill="white", colour="white"), strip.text=element_text(face=NULL, size=7), axis.text=element_text(size=5), 
        legend.position="bottom", legend.title=element_blank())
ggsave("figures/statin_comp_acc_treatgroups_4b.pdf", width=9, height=6.25, scale=1, units="cm")

# for nocovars version:
my.colors <- c("royalblue1", "seagreen3")
# ...
ggsave("figures/statin_comp_acc_nocovars_treatgroups_supp3a.pdf", width=12, height=7, scale=1, units="cm")

