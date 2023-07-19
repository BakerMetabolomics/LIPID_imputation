
##############################################################################
##
##  Validate with multi imputed MEASURED LIPID species: CVD-lipid associations
##
##  Association analysis to be done at three levels:
##  1. original measured lipids in LIPID
##  2. imputed measured lipids (exact)
##  3. multiple imputed measured lipids (average association)
##
##  Corresponds to results in Supplementary Figure 3.
##############################################################################

library("tidyverse")
library("magrittr")
library("broom")
library("openxlsx")
library("RColorBrewer")

##############################################################################
##
##  Load saved multiple imputed MEASURED only LIPID species
##
##############################################################################

multi_LIPID_predict <- readRDS('results/multi_LIPID_predicted_measured_species_294.rds')

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



##############################################################################
##############################################################################
##
##  Select original and multiple imputed MEASURED lipid species for analysis
##
##############################################################################
##############################################################################



##############################################################################
##  Prepare LIPID trial data
##############################################################################

## Note: few lipid names bellow will need to be changed to be in line with matching AusDiab names from 12/2021
lipid <- read_csv(file="/Volumes/labs/Metabolomics/Projects/LIPID/Databases/2017_09_12 - Lipid Trial Database - Mapped to Ausdiab.csv", 
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf) %>% 
  rename(`PC(39:5) (b)`=`PC (39:5) (b)`, `PE(O-18:0/22:5)`=`PE(O-18:0/22:5) (a)`, `PC(28:0)`=`PC 28:0`)

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
    select(id, eventCVD), 
  by="id") %>% 
  relocate(eventCVD, .after=treat)

## select only followup, exclude 13 discordant lipids and add cvd variable
lipid_composite_flup_cvd <- left_join(
  lipid_composite %>% 
    filter(id %in% followup_id) %>% 
    select(-all_of(exclude_pcor_13)), 
  lipid %>% 
    mutate(id=as.character(ID)) %>% 
    filter(id %in% followup_id) %>% 
    select(id, eventCVD), 
  by="id") %>% 
  relocate(eventCVD, .after=treat)



##############################################################################
##############################################################################
##
##  Univariate association of original & imputed lipids with CVD outcome 
##
##############################################################################
##############################################################################


##############################################################################
## Function for checking various univariate associations
##############################################################################

## Odds ratios, CI results
TestUnivarAssociations_CatOutcome_Exp <- function(data.in, outcome) {
  
  fo <- paste(outcome, "~")
  
  result <- data.in %>% 
    select(-c(id, all_of(outcome))) %>% 
    map(.f= ~ glm(formula=as.formula(paste(fo, ".x")), family="binomial", data=data.in)) %>% 
    map(.f=tidy, conf.int=T, conf.level=0.95, exponentiate=T) %>% 
    map(.f=select, term, estimate, conf.low, conf.high) %>% 
    map_dfr(.f=slice_tail, n=1, .id="lipid") %>% 
    select(lipid, estimate, conf.low, conf.high) %>% 
    rename(OR=estimate)
  
  return(result)
  
} 

## Betas, p value results
TestUnivarAssociations_CatOutcome_Orig <- function(data.in, outcome) {
  
  fo <- paste(outcome, "~")
  
  result <- data.in %>% 
    select(-c(id, all_of(outcome))) %>% 
    map(.f= ~ glm(formula=as.formula(paste(fo, ".x")), family="binomial", data=data.in)) %>% 
    map(.f=tidy, conf.int=T, conf.level=0.95, exponentiate=F) %>% 
    map(.f=select, term, estimate, std.error, p.value) %>% 
    map_dfr(.f=slice_tail, n=1, .id="lipid") %>% 
    select(lipid, estimate, std.error, p.value)
  
  return(result)
  
} 


TestUnivarAssociations_CatOutcome <- function(data.in, outcome) {
  
  or.res <- TestUnivarAssociations_CatOutcome_Exp(data.in, outcome)
  beta.res <- TestUnivarAssociations_CatOutcome_Orig(data.in, outcome)
  res <- left_join(or.res, beta.res)
  
}


###############################################################################
## Function to do the job on 10 x 1-3 data versions 
###############################################################################
GlmMultiEstimateCat <- function(orig_data, multi_imputed_data, outcome_var) {
  
  ## set variable selectors
  lipid_names <- orig_data %>% select(`CE(14:0)`:`TG(58:8) [NL-22:6]`) %>% colnames()
  
  original <- orig_data %>% 
    select(id, all_of(outcome_var), `CE(14:0)`:`TG(58:8) [NL-22:6]`)
  
  ## set container
  result_3x_10 <- list()
  
  ## loop
  for(j_ in seq(length(multi_imputed_data))) {
    
    # set imputed data only
    imputed <- left_join(original %>% 
                           select(id, all_of(outcome_var)), 
                         multi_imputed_data[[j_]]$LIPID_miss %>% 
                           mutate(across(.cols=`CE(14:0)`:`TG(58:8) [NL-22:6]`, .fns=scale)), 
                         by="id")
    
    # set stochastic imputed data only
    imputed_sto <- left_join(original %>% 
                               select(id, all_of(outcome_var)), 
                             multi_imputed_data[[j_]]$LIPID_miss_sto %>% 
                               mutate(across(.cols=`CE(14:0)`:`TG(58:8) [NL-22:6]`, .fns=scale)), 
                             by="id")
    
    if (j_ == 1) {
      
      # bind all in a list
      data_3x <- list('original'=original,
                      'imputed'=imputed,
                      'imputed_sto'=imputed_sto)
      
      # do the analysis
      result <- map(.x=data_3x, .f=TestUnivarAssociations_CatOutcome, outcome=outcome_var)
      
    } else {
      
      # bind all in a list
      data_1x <- list('imputed_sto'=imputed_sto)
      
      
      # do the analysis
      result <- map(.x=data_1x, .f=TestUnivarAssociations_CatOutcome, outcome=outcome_var)
      
    }
    
    result_3x_10[[j_]] <- result
    
  }
  
  return(result_3x_10)
  
}



###############################################################################
## Estimate CVD - lipid associations for original, imputed, & 10x sto. imputed
## with BASELINE lipids
###############################################################################
cvd_lipid_or_result_base <- GlmMultiEstimateCat(orig_data=lipid_composite_base_cvd, 
                                           multi_imputed_data=multi_LIPID_predict, 
                                           outcome_var="eventCVD")

saveRDS(cvd_lipid_or_result_base, 'results/validation_multi_impute_cvd_measured_lipid_OR_results_base.rds')
# cvd_lipid_or_result_base <- readRDS('results/validation_multi_impute_cvd_measured_lipid_OR_results_base.rds')


###############################################################################
## Summarise the results

# original and imputed lipid associations are already summarised
# further pooled summary required for the multiple imputations

## Original lipids
orig_result_summ_base <- cvd_lipid_or_result_base[[1]]$original %>% 
  rename(OR_low=conf.low, OR_high=conf.high, se=std.error, p=p.value) %>% 
  mutate(p_BH=p.adjust(p, method="BH"))

## Exact imputed lipids
impu_result_summ_base <- cvd_lipid_or_result_base[[1]]$imputed %>% 
  rename(OR_low=conf.low, OR_high=conf.high, se=std.error, p=p.value) %>% 
  mutate(p_BH=p.adjust(p, method="BH"))

## Multiple imputed lipids; Rubin's rules for pooling estimates
mi_result_base <- bind_rows(map(cvd_lipid_or_result_base, 'imputed_sto')) %>% 
  arrange(lipid)

mi_result_summ_base <- mi_result_base %>% 
  group_by(lipid) %>% 
  mutate(m_estimate=mean(estimate)) %>% 
  mutate(var_w=mean(std.error^2)) %>% 
  mutate(var_b=sum((estimate-m_estimate)^2)/9) %>% 
  mutate(m_var=var_w+var_w*1.1) %>% 
  mutate(m_se=sqrt(m_var)) %>% 
  mutate(m_est_lo=m_estimate-qnorm(p=0.975)*m_se) %>% 
  mutate(m_est_hi=m_estimate+qnorm(p=0.975)*m_se) %>% 
  mutate(m_p=2*pnorm(-abs(m_estimate/m_se))) %>% 
  mutate(m_OR=exp(m_estimate)) %>% 
  mutate(m_OR_low=exp(m_est_lo)) %>% 
  mutate(m_OR_high=exp(m_est_hi)) %>% 
  summarise(OR=mean(m_OR), OR_low=mean(m_OR_low), OR_high=mean(m_OR_high), estimate=mean(m_estimate), se=mean(m_se), p=mean(m_p)) %>% 
  mutate(p_BH=p.adjust(p, method="BH"))


## Put all summaries together; Create lipid classes
cvd_lipid_assoc_summary_base <- bind_rows(list("measured"=orig_result_summ_base, "imputed"=impu_result_summ_base, "multiple_imputed"=mi_result_summ_base), 
                                     .id="data")

cvd_lipid_assoc_summary_base <- cvd_lipid_assoc_summary_base %>% 
  mutate(class=case_when(
    str_starts(lipid, "PC\\(O") ~ "PC(O)", 
    str_starts(lipid, "PC\\(P") ~ "PC(P)", 
    str_starts(lipid, "PE\\(O") ~ "PE(O)", 
    str_starts(lipid, "PE\\(P") ~ "PE(P)", 
    str_starts(lipid, "LPC\\(O") ~ "LPC(O)", 
    str_starts(lipid, "\\w+\\(") ~ str_extract(lipid, "\\w+"), 
    str_starts(lipid, "\\w+") ~ str_extract(lipid, "\\w+")
  )) %>% 
  relocate(class, .after=lipid)


## List of validation summary tables ----------------------------------------##
cvd_lipid_assoc_summ_tables_base <- list('measured'=cvd_lipid_assoc_summary_base %>% filter(data=="measured"), 
                                    'imputed'=cvd_lipid_assoc_summary_base %>% filter(data=="imputed"), 
                                    'multiple imputed'=cvd_lipid_assoc_summary_base %>% filter(data=="multiple_imputed"), 
                                    'all together'=cvd_lipid_assoc_summary_base)

## Write a multi-sheet excel file
write.xlsx(cvd_lipid_assoc_summ_tables_base, file="results/validation_multi_impute_cvd_measured_lipid_OR_summary_base.xlsx")

## Create plots in Excel (Thy)



