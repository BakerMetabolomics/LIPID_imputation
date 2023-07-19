
##############################################################################
##
## Find irreconcilable species between AusDb & LIPID to remove from predictors 
##
## Corresponds to results in Figure 1.
##############################################################################

library("tidyverse")
library("magrittr")
library("pdist")
library("corpcor")

##############################################################################
##
##  Load saved working datasets 
##
##############################################################################

## Load saved composite datasets
ausdb_composite <- readRDS("data_derived/Ausdiab_composite_working_data_log_scale.rds") # reference data frame. variables: id, age, sex, bmi, chol_treat, 829 lipid species
lipid_composite <- readRDS("data_derived/LIPID_composite_working_data_log_scale.rds") # target data frame. variables: id, age, sex, bmi, treat_12, 307 lipid species

## Load splits of composite analysis lipid names 
Ausdiab_LIPID_lipid_names_collection <- readRDS('data_derived/Ausdiab_LIPID_lipid_names_collection.rds')

## Select matching variables needed here
ausdb307 <- ausdb_composite %>% select(id, all_of(Ausdiab_LIPID_lipid_names_collection$matching_307c))
lipid307 <- lipid_composite %>% select(id, all_of(Ausdiab_LIPID_lipid_names_collection$matching_307c))

## Make matrices for glmnet
ausdb307_mat <- as.matrix(ausdb307[,-1])
lipid307_mat <- as.matrix(lipid307[,-1])



###############################################################################
##
## Remove discrepant lipids between two datasets one-by-one based on differences 
## in their Partial correlation vectors. Run pcor after every removal, remove
## the most discrepant lipid, stop when there's no apparent outliers (3sd, 3mad)
##
###############################################################################


## Make a function for this task --------------------------------------------##

RemoveDissimilarLipidsPcor <- function(mat1, mat2, FUN, num) {
  
  # calculate partial correlations in two data sets
  pcor1 <- pcor.shrink(mat1, lambda=0, verbose=F)
  diag(pcor1) <- 0
  
  pcor2 <- pcor.shrink(mat2, lambda=0, verbose=F) 
  diag(pcor2) <- 0
  
  # calculate pcor distance between corresponding lipids in two datasets
  pcor_dist <- pdist(pcor1, pcor2) %>% as.matrix() %>% diag()
  names(pcor_dist) <- colnames(pcor1)
  
  # find how many dissimilar lipids there are (> x sd, mad variation)
  med_x_var <- median(pcor_dist) + num*match.fun(FUN)(pcor_dist)
  n_dissimilar <- length(pcor_dist[pcor_dist > med_x_var])
  
  fro_ave <- sum(pcor_dist^2) / (length(pcor_dist)^2)
  
  print( paste0("Initial # dissim.=", n_dissimilar, " Ave. Frobenius=", fro_ave ))
  
  # containers
  excluded <- vector("character")
  pcor_dist_series <- list(pcor_dist)
  frobenius_ave <- c(fro_ave)
  count <- 1
  
  repeat {
    
    worst <- which(pcor_dist == max(pcor_dist))
    worst_name <- names(worst)
    print(paste0("Excluding ", count, ":  ", worst_name))
    
    mat1 <- mat1[, -worst]
    mat2 <- mat2[, -worst]
    
    # calculate partial correlations in reduced data sets
    pcor1 <- pcor.shrink(mat1, lambda=0, verbose=F)
    diag(pcor1) <- 0
    
    pcor2 <- pcor.shrink(mat2, lambda=0, verbose=F) 
    diag(pcor2) <- 0
    
    # calculate new pcor distance between corresponding lipids
    pcor_dist <- pdist(pcor1, pcor2) %>% as.matrix() %>% diag() # Euclid dist. (sqrt of summed squared dists.)
    names(pcor_dist) <- colnames(pcor1)
    
    # find how many dissimilar lipids there are now (> x sd, mad variation)
    med_x_var <- median(pcor_dist) + num*match.fun(FUN)(pcor_dist)
    n_dissimilar <- length(pcor_dist[pcor_dist > med_x_var])
    
    fro_ave <- sum(pcor_dist^2) / (length(pcor_dist)^2) # mean squared Frobenius dist.
    
    print( paste0("# dissim.=", n_dissimilar, " Ave. Frobenius: ", fro_ave ))
    
    excluded <- append(excluded, worst_name)
    pcor_dist_series <- append(pcor_dist_series, list(pcor_dist))
    frobenius_ave <- append(frobenius_ave, fro_ave)
    count <- count + 1
    
    if (n_dissimilar == 0) break
    
  }
  
  remaining <- colnames(mat1)
  return(list("excluded"=excluded, "remaining"=remaining, "pcor_dist_series"=pcor_dist_series, "frobenius_ave"=frobenius_ave))
  
}


## Run the function with few variations
lipids_to_exclude_3mad <- RemoveDissimilarLipidsPcor(ausdb307_mat, lipid307_mat, mad, 3) # 9
lipids_to_exclude_2.5mad <- RemoveDissimilarLipidsPcor(ausdb307_mat, lipid307_mat, mad, 2.5) # 13
lipids_to_exclude_2mad <- RemoveDissimilarLipidsPcor(ausdb307_mat, lipid307_mat, mad, 2) # 21


exclude_predictors_pcor_measure <- list('mad_2'=lipids_to_exclude_2mad, 'sd_2'=lipids_to_exclude_2sd)
saveRDS(exclude_predictors_pcor_measure, 'results/exclude_predictors_pcor_measure.rds')
# exclude_predictors_pcor_measure <- readRDS('results/exclude_predictors_pcor_measure.rds')



#-----------------------------------------------------------------------------#
## Plots

## Plot overall Frobenius distances -----------------------------------------##
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
ggplot(tibble(x=0:21, y=exclude_predictors_pcor_measure$mad_2$frobenius_ave), aes(x, y)) +
  geom_line(size=0.4) + 
  geom_point(size=1.75, shape=21, fill="tomato", colour="white", stroke=0.75) +
  theme_bw(base_size=8) + 
  scale_y_continuous(label=scientific_10) + 
  labs(x="Number of excluded lipids", y="Average squared distance") + 
  geom_vline(aes(xintercept=9), linetype="dotted", size=0.4) + 
  geom_vline(aes(xintercept=13), linetype="dotted", size=0.4) + 
  theme(axis.text=element_text(size=6))   
ggsave("figures/adj_average_pcor_sq_dist.pdf", width=8, height=6, scale=1, units="cm")
ggsave("figures/adj_average_pcor_sq_dist.tiff", width=8, height=6, scale=1, units="cm")
ggsave("figures/adj_average_pcor_sq_dist.svg", width=8, height=6, scale=1, units="cm")

scales::scientific
scale_y_continuous(label=function(x) format(x, scientific = TRUE))

## Plot distribution of per lipid distances for the first 15 exclusions (lipids_to_exclude_2.5mad)
dist_df <- tibble(dist=unlist(exclude_predictors_pcor_measure$mad_2$pcor_dist_series[1:16])^2, 
                  excluded=rep(0:15, times=307:(307-15))) %>% 
  mutate(excluded=factor(as.character(excluded), levels=as.character(0:15)))

ggplot(dist_df, aes(dist)) +
  geom_histogram(binwidth=0.03, colour="white", size=0.003, fill="royalblue3") + 
  theme_bw(base_size=7.5) + 
  labs(x="Squared distance") + 
  scale_y_continuous(breaks=c(0,25,50)) + 
  facet_wrap(~excluded, scales="fixed", ncol=4) + 
  theme(axis.text=element_text(size=5)) 
ggsave("figures/adj_distribution_of_pcor_sq_dist.pdf", width=10, height=8, scale=1, units="cm")
ggsave("figures/adj_distribution_of_pcor_sq_dist.tiff", width=10, height=8, scale=1, units="cm")
ggsave("figures/adj_distribution_of_pcor_sq_dist.svg", width=10, height=8, scale=1, units="cm")


## Table excluded predictors
excluded_predictors <- tibble(lipids_to_exclude_pcor=exclude_predictors_pcor_measure$mad_2$excluded[1:13])
write.xlsx(excluded_predictors, file = "results/excluded_predictors_table.xlsx")


#-----------------------------------------------------------------------------#
# Some checks:

# calculate partial correlations in two data sets before any removals
pcor_ausdb <- pcor.shrink(ausdb307_mat, lambda=0, verbose=F)
diag(pcor_ausdb) <- 0

pcor_lipid <- pcor.shrink(lipid307_mat, lambda=0, verbose=F) 
diag(pcor_lipid) <- 0

# Check some individual lipid species
xbn <- "PC(O-36:3)"
xb <- which(colnames(ausdb307_mat)==xbn)

# pcor plot ausdb vs. lipid
plot(pcor_ausdb[,xbn], pcor_lipid[,xbn], xlim=c(-0.5,1), ylim=c(-0.5,1), main=xbn)

# closest in ausdb
ybn <- names(which.max(pcor_ausdb[,xbn]))

plot(ausdb307_mat[,xbn], ausdb307_mat[,ybn])
plot(lipid307_mat[,xbn], lipid307_mat[,ybn])

# closest in lipid
zbn <- names(which.max(pcor_lipid[,xbn]))

plot(ausdb307_mat[,xbn], ausdb307_mat[,zbn])
plot(lipid307_mat[,xbn], lipid307_mat[,zbn])


