
##############################################################################
##
##  Make working version of AusDiab & LIPID inc. composite species definitions
##
##  Corresponds to results in Supplementary Tables 1-3.
##############################################################################

library("tidyverse")
library("magrittr")
library("readxl")
library("openxlsx")
library("tableone")


##############################################################################
##
##  Load AusDiab and prepare variables
##
##############################################################################

ausdb <- read_csv(file="/Volumes/labs/Metabolomics/Projects/Ausdiab/Databases/2021_12_16 - Ausdiab lipidomics and clinical data combined.csv", 
                  col_select = c(id, drage_00, drsex_00_n, bmi_00, drfatper_00, systolic_00, diastoli_00, q23_tabl_00, 
                                 q20_angi_00, q20_coro_00, q20_stro_00, CVDdnfev_10, cvd_d, cva_d, choltabl_00, 
                                 `AC(12:0)`:`Ubiquinone`), 
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf)


# Some fixes: select few extra variables, fix variable names, delete redundant ones
ausdb <- ausdb %>% 
  mutate(chol_treat=if_else(choltabl_00=="Yes", 1, 0, missing=0)) %>% 
  mutate(bp_treat=if_else(q23_tabl_00=="yes", 1, 0, missing=0)) %>% 
  mutate(base_cvd=if_else(q20_angi_00=="yes" | q20_coro_00=="yes" | q20_stro_00=="yes", 1, 0, missing=0)) %>% 
  mutate(inci_cvd=if_else(CVDdnfev_10=="CVD" | cvd_d=="CVD death on or before 30Nov2013" | cva_d=="CVA death on or before 30Nov2013", 1, 0, missing=0)) %>% 
  mutate(comb_cvd=if_else(base_cvd==1 | inci_cvd==1, 1, 0)) %>% 
  rename(age=drage_00, bmi=bmi_00, sex=drsex_00_n) %>% 
  mutate(id=as.character(id)) %>% 
  select(id, age, sex, bmi, chol_treat, bp_treat, base_cvd, inci_cvd, comb_cvd, `AC(12:0)`:`Ubiquinone`)



##############################################################################
##
##  Load individual & composite lipids mapping files 
##
##############################################################################

## Individual lipids mapping.
LIPID_map <- read_xlsx("/Volumes/labs/Metabolomics/Projects/LIPID/Analysis/2021_07_07 - LIPID mapping NR/2022_01_20 - LIPID mapped to Ausdiab_surrogates_updates.xlsx", 
                       sheet=1, range ="A1:D343", na=c("", "#N/A"))

## Composite lipids mapping.
composite_mapping <- read_csv("/Volumes/labs/Metabolomics/Projects/LIPID/Analysis/2021_07_07 - LIPID mapping NR/2022_03_29 - Sum composition updated KH.csv", 
                              trim_ws=TRUE) 


## Get AusDb lipid names
ausdb_lipid_names <- colnames(ausdb)[10:756]

## Create individual matching/missing lipid categories (names)
matching_lipid_names <- ausdb_lipid_names[ausdb_lipid_names %in% LIPID_map$`LIPID possible names - Thy's match`] # 307 matches
missing_lipid_names <- setdiff(ausdb_lipid_names, matching_lipid_names) # 440 to impute

## Check what isn't mapped
nomap_lipid_names <- LIPID_map %>% 
  filter(!`LIPID possible names - Thy's match` %in% ausdb_lipid_names) %>% 
  pull(`LIPID possible names - Thy's match`) # 35, mostly some DG, PC, TG species

## Get lipid names used in creation of composite lipids
noncomposite_names <- vector("character")
for (i in seq(nrow(composite_mapping))) {
  noncomposite_names <- append(noncomposite_names, str_split(composite_mapping$Ausdiab[i], pattern=", ")[[1]])
} # 176 lipids

noncomposite_matching_names <- intersect(noncomposite_names, matching_lipid_names) # 82
noncomposite_missing_names <- intersect(noncomposite_names, missing_lipid_names) # 94

## Some checks
sum(noncomposite_names %in% colnames(ausdb)) # 176 = 82 + 94



##############################################################################
##
##  Load LIPID and prepare variables / timepoints
##
##############################################################################

lipid <- read_csv(file="/Volumes/labs/Metabolomics/Projects/LIPID/Databases/2017_09_12 - Lipid Trial Database - Mapped to Ausdiab.csv", 
                  col_select = c(ID, pid, period, TREAT, age, sexn, bmi, `CE(14:0)`:`Hex3Cer(d18:1/24:1)`), 
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf)

## Some fixes: adjust treatment variable, remove period, pid for now
lipid <- lipid %>% 
  rename(sex=sexn, id=ID) %>% 
  mutate(treat_12=if_else(TREAT==1 & period==12, 1, 0)) %>% 
  mutate(id=as.character(id)) %>% 
  select(-c(TREAT, period, pid)) # remove period, pid as required

## Fix few lipid names, spaces, and select only matching species (for the purpose of prediction)
lipid <- lipid %>% 
  rename(`PC(39:5) (b)`=`PC (39:5) (b)`, `PE(O-18:0/22:5)`=`PE(O-18:0/22:5) (a)`, `PC(28:0)`=`PC 28:0`) %>% 
  select(id, age, sex, bmi, treat_12, all_of(matching_lipid_names))

## Some checks
sum(noncomposite_names %in% colnames(lipid)) ## 82
all.equal(intersect(noncomposite_names, colnames(lipid)), intersect(noncomposite_names, matching_lipid_names)) # TRUE

## in both data sets 0=female, 1=male



##############################################################################
##
##  Change to composite LIPID lipid names to match AusDb composite names
##
##############################################################################

## Make composite / noncomposite LIPID matches
LIPID_matches <- tibble(composite=vector("character"), noncomposite=vector("character"))
for (i in seq(nrow(composite_mapping))) {
  LIPID_matches[i, 1] <- composite_mapping$Sum_comp[i]
  LIPID_matches[i, 2] <- intersect(colnames(lipid)[6:312], str_split(composite_mapping$Ausdiab[i], pattern=", ")[[1]])
}

## Change LIPID lipid names (back) to composite
LIPID_lnames <- tibble(ind=colnames(lipid), comp=colnames(lipid))
for (i in 1:312) {
  if (LIPID_lnames$ind[i] %in% LIPID_matches$noncomposite) {
    j <- which(LIPID_matches$noncomposite %in% LIPID_lnames$ind[i])
    LIPID_lnames$comp[i] <- LIPID_matches$composite[j]
  }
}

colnames(lipid) <- LIPID_lnames$comp

## NOTE: 
## At point of imputed data, names will go back to non-composite to match AusDiab



##############################################################################
##############################################################################
##
##  Generate composite Ausdiab lipids
##
##############################################################################
##############################################################################

## Sum up original lipids according to composite_mapping file
ausdb_composite <- ausdb
for (i in seq(nrow(composite_mapping))) {
  
  comp_name <- composite_mapping$Sum_comp[i]
  ind_names <- str_split(composite_mapping$Ausdiab[i], pattern=", ")[[1]]
  #i_names <- which(colnames(dat) %in% ind_names)
  
  if (length(ind_names)==1) {
    
    ausdb_composite <- ausdb_composite %>% 
      mutate("{comp_name}" := .data[[ind_names[1]]], .keep="all")
    
  } else if (length(ind_names)==2) {
    
    ausdb_composite <- ausdb_composite %>% 
      mutate("{comp_name}" := .data[[ind_names[1]]]+.data[[ind_names[2]]], .keep="all")
    
  } else if (length(ind_names)==3) {
    
    ausdb_composite <- ausdb_composite %>% 
      mutate("{comp_name}" := .data[[ind_names[1]]]+.data[[ind_names[2]]]+.data[[ind_names[3]]], .keep="all")
    
  } else if (length(ind_names)==4) {
    
    ausdb_composite <- ausdb_composite %>% 
      mutate("{comp_name}" := .data[[ind_names[1]]]+.data[[ind_names[2]]]+.data[[ind_names[3]]]+.data[[ind_names[4]]], .keep="all")
    
  } 
  
}


## Create composite matching/missing lipid categories (names)
## (keeping all non-composite lipids in "missing" category to be predicted)
missing_522c <- union(missing_lipid_names, noncomposite_matching_names) %>% sort()
matching_307c <- setdiff(colnames(ausdb_composite)[10:838], missing_522c) %>% sort()


## Save a handy collection of lipid name categories
Ausdiab_LIPID_lipid_names_collection <- list(matching_307i=matching_lipid_names, missing_440i=missing_lipid_names, 
                                             matching_307c=matching_307c, missing_522c=missing_522c, 
                                             noncomposite_matching_82=noncomposite_matching_names, noncomposite_missing_94=noncomposite_missing_names, 
                                             nomap_lipid_names_35=nomap_lipid_names)
saveRDS(Ausdiab_LIPID_lipid_names_collection, 'data_derived/Ausdiab_LIPID_lipid_names_collection.rds')
# Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')



##############################################################################
##
##  Transform, save reduced datasets for the imputation workflow
##
##############################################################################

## Log transform and scale (regardless of glmnet, required for test set compatibility)
ausdb_composite <- ausdb_composite %>% 
  mutate(across(.cols=c(bmi, `AC(12:0)`:`SM(41:2)`), .fns=log)) %>% 
  mutate(across(.cols=`AC(12:0)`:`SM(41:2)`, .fns=scale))

lipid_composite <- lipid %>% 
  mutate(across(.cols=c(bmi, `CE(14:0)`:`TG(58:8) [NL-22:6]`), .fns=log)) %>% 
  mutate(across(.cols=`CE(14:0)`:`TG(58:8) [NL-22:6]`, .fns=scale))

saveRDS(ausdb_composite, 'data_derived/Ausdiab_composite_working_data_log_scale.rds') # 9922 x 838
# ausdb_composite <- readRDS("data_derived/Ausdiab_composite_working_data_log_scale.rds")
saveRDS(lipid_composite, 'data_derived/LIPID_composite_working_data_log_scale.rds') # 11773 x 312
# lipid_composite <- readRDS("data_derived/LIPID_composite_working_data_log_scale.rds")


rm(ausdb, lipid)



##############################################################################
##
##  Table composite mapping
##
##############################################################################

comp_map <- composite_mapping %>% 
  select(Sum_comp, Ausdiab) %>% 
  mutate("Number of summed lipid species"=str_count(Ausdiab, ",")+1) %>% 
  rename("Composite lipid"=Sum_comp, "Summed lipid species (AusDiab)"=Ausdiab) %>% 
  select("Composite lipid", "Number of summed lipid species", "Summed lipid species (AusDiab)")

write.xlsx(comp_map, file = "results/composite_mapping_table.xlsx")

comp_map %>% filter(`Number of summed lipid species`==2) # 70
comp_map %>% filter(`Number of summed lipid species`==3) # 12



##############################################################################
##
##  Table 1. Anthropometric & clinical characteristics, AusDiab & LIPID
##
##############################################################################

## Starting from the whole data sets
ausdb <- read_csv(file="/Volumes/labs/Metabolomics/Projects/Ausdiab/Databases/2021_12_16 - Ausdiab lipidomics and clinical data combined.csv",  
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf)

lipid <- read_csv(file="/Volumes/labs/Metabolomics/Projects/LIPID/Databases/2017_09_12 - Lipid Trial Database - Mapped to Ausdiab.csv", 
                  na = c("", "NA"), trim_ws = TRUE, guess_max = Inf)

## Focus on age, sex, BMI, chol, HDL, LDL, TG, lipid lowering med., diabetes, CVD preval/incid
ausdb <- ausdb %>% 
  select(drage_00, drsex_00, bmi_00, chol_00, hdl_00, ldl_00, trig_00, choltabl_00, diabstat_00, q20_angi_00, q20_coro_00, q20_stro_00, CVDdnfev_10) %>% 
  mutate("lipid treatment"=if_else(choltabl_00=="Yes", "Yes", "No", missing="No")) %>% 
  mutate("diabetes"=if_else(diabstat_00=="kdm" | diabstat_00=="New DM", "Yes", "No")) %>% 
  mutate("previous CVD"=if_else(q20_angi_00=="yes" | q20_coro_00=="yes" | q20_stro_00=="yes", "Yes", "No", missing="No")) %>% 
  mutate("incident CVD"=if_else(CVDdnfev_10=="CVD", "Yes", "No", missing="No")) %>% 
  rename(age=drage_00, gender=drsex_00, BMI=bmi_00, "total cholesterol"=chol_00, HDL=hdl_00, LDL=ldl_00, triglycerides=trig_00) %>% 
  select(age, gender, BMI, "total cholesterol", HDL, LDL, triglycerides, "lipid treatment", diabetes, "previous CVD", "incident CVD")
  
lipid <- lipid %>% 
  select(period, age, sexn, bmi, chol0q, hdl0q, ldl0q, trig0q, TREAT, diab, eventCVD) %>% 
  filter(period==0) %>% 
  mutate(gender=if_else(sexn==1, "Male", "Female")) %>% 
  mutate("lipid treatment"=if_else(TREAT==1, "Yes", "No")) %>% 
  mutate("diabetes"=if_else(diab==1, "Yes", "No")) %>% 
  mutate("previous CVD"="Yes") %>% 
  mutate("incident CVD"=if_else(eventCVD==1, "Yes", "No")) %>% 
  rename(BMI=bmi, "total cholesterol"=chol0q, HDL=hdl0q, LDL=ldl0q, triglycerides=trig0q) %>% 
  select(age, gender, BMI, "total cholesterol", HDL, LDL, triglycerides, "lipid treatment", diabetes, "previous CVD", "incident CVD")
  
## put them together for table
joint_table <- bind_rows(list("AusDiab"=ausdb, "LIPID"=lipid), .id="study") %>% 
  mutate(across(where(is.character), factor))


## TableOne version #######################################
table_one <- CreateTableOne(vars=c("age", "gender", "BMI", "total cholesterol", "HDL", "LDL", "triglycerides", "lipid treatment", "diabetes", "previous CVD", "incident CVD"), 
                            strata="study",
                            data=joint_table,
                            test=FALSE)
# 1 lvl shown
table_one_mat <- print(table_one, 
                       explain=TRUE,
                       printToggle=TRUE, 
                       test=FALSE,
                       noSpaces=TRUE,
                       showAllLevels=FALSE)
write.csv(table_one_mat, "results/anthropo_tableone_1lvl.csv")




###########################################################
## Creating Tables for publication: lipid names breakdown
###########################################################

##-------------------------------------------------------##
## List for lipid names
LIPID_AusDb_names <- list("LIPID_1-to-1_match"=setdiff(Ausdiab_LIPID_lipid_names_collection$matching_307i, 
                                                       Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82), 
                          "LIPID_1-to-23_match"=Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82, 
                          "LIPID_no_match"=Ausdiab_LIPID_lipid_names_collection$nomap_lipid_names_35, 
                          "AusDb_1-to-1_match"=setdiff(Ausdiab_LIPID_lipid_names_collection$matching_307i, 
                                                       Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82), 
                          "AusDb_23-to-1_match"=c(Ausdiab_LIPID_lipid_names_collection$noncomposite_matching_82, 
                                                   Ausdiab_LIPID_lipid_names_collection$noncomposite_missing_94), 
                          "AusDb_no_match"=setdiff(Ausdiab_LIPID_lipid_names_collection$missing_440i, 
                                                   Ausdiab_LIPID_lipid_names_collection$noncomposite_missing_94))

## Write a multi-sheet excel file
write.xlsx(LIPID_AusDb_names, file="results/LIPID_AusDb_lipid_name_categories.xlsx")


##-------------------------------------------------------##
## List for lipid mappings
LIPID_map_categories <- list("LIPID_1-to-1_match"=LIPID_map %>% 
                               filter(`LIPID possible names - Thy's match` %in% LIPID_AusDb_names$`LIPID_1-to-1_match`), 
                             "LIPID_1-to-23_match"=LIPID_map %>% 
                               filter(`LIPID possible names - Thy's match` %in% LIPID_AusDb_names$`LIPID_1-to-23_match`), 
                             "LIPID_no_match"=LIPID_map %>% 
                               filter(`LIPID possible names - Thy's match` %in% LIPID_AusDb_names$`LIPID_no_match`))

## Write a multi-sheet excel file
write.xlsx(LIPID_map_categories, file="results/LIPID_map_lipid_name_categories.xlsx")


##-------------------------------------------------------##
## Another LIPID_to_AusDiab_map view (from manually created Excel sheet)
LIPID_to_AusDiab_map <- read_xlsx("data_derived/LIPID_to_AusDiab_map_2022_12_12_AD.xlsx", 
                       sheet=1, range ="A1:E343", na=c("", "#N/A"))

write.csv(LIPID_to_AusDiab_map, "data_derived/LIPID_to_AusDiab_map_2022_12_12_AD.csv")

# it was then manually checked by Nat


##-------------------------------------------------------##
## Another AusDiab_to_AusDiab_composite name map view 
ausdb_to_comp_ausdb <- composite_mapping %>% 
  separate_rows(Ausdiab, sep=", ") %>% 
  select(Ausdiab, Sum_comp) %>% 
  group_by(Sum_comp) %>% mutate(n=n())

ausdb_to_ausdb <- tibble(Ausdiab=ausdb_lipid_names, Sum_comp=ausdb_lipid_names, n=1) %>% 
  filter(!Ausdiab %in% ausdb_to_comp_ausdb$Ausdiab) %>% 
  mutate(n=if_else(!Sum_comp %in% Ausdiab_LIPID_lipid_names_collection$matching_307c, 0, 1))

ausdb_map <- bind_rows(ausdb_to_ausdb, ausdb_to_comp_ausdb) %>% 
  arrange(Sum_comp) %>% 
  mutate(Grouping=case_when(n==0~"No LIPID match", n==1~"Single", n>=2~"Composite")) %>% 
  select(Ausdiab, Sum_comp, Grouping, n) %>% 
  rename(AusDiab=Ausdiab, "Single or Composite interim name"=Sum_comp, "Number of matching AusDiab species"=n)

write.csv(ausdb_map, "data_derived/AusDiab_map_2022_12_12_AD.csv", row.names=FALSE)
write.xlsx(ausdb_map, file="data_derived/AusDiab_map_2022_12_12_AD.xlsx")

