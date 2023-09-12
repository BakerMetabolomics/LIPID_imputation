
##############################################################################
##
##  Simulate data to evaluate code
##
##  Output names correspond to files used in scripts 1-6.
##############################################################################

library(openxlsx)
library(MASS)
library(Matrix)

# Set the seed for reproducibility
set.seed(6495)

# Specify the number of subjects to simulate
n <- 1000

##############################################################################
##
##  Create Ausdiab Dataset
##
##############################################################################

# Create a data frame with completely simulated data
data_Ausdiab <- data.frame(
    id = seq(1, n),
    drage_00 = sample(20:80, n, replace = TRUE),
    drsex_00_n = sample(c("Male", "Female"), n, replace = TRUE),
    bmi_00 = rnorm(n, mean = 25, sd = 5),
    drfatper_00 = rnorm(n, mean = 25, sd = 5),
    systolic_00 = rnorm(n, mean = 120, sd = 10),
    diastoli_00 = rnorm(n, mean = 80, sd = 10),
    q23_tabl_00 = sample(c("yes", "no"), n, replace = TRUE, prob = c(0.2, 0.8)),
    q20_angi_00 = sample(c("yes", "no"), n, replace = TRUE), prob = c(0.2, 0.8),
    q20_coro_00 = sample(c("yes", "no"), n, replace = TRUE, prob = c(0.2, 0.8)),
    q20_stro_00 = sample(c("yes", "no"), n, replace = TRUE, prob = c(0.2, 0.8)),
    CVDdnfev_10 = sample(c("CVD", "No CVD"), n, replace = TRUE, prob = c(0.2, 0.8)),
    cvd_d = sample(c("CVD death on or before 30Nov2013", "No CVD death"), n, replace = TRUE, prob = c(0.2, 0.8)),
    cva_d = sample(c("CVA death on or before 30Nov2013", "No CVA death"), n, replace = TRUE, prob = c(0.2, 0.8)),
    choltabl_00 = sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.2, 0.8))
)

# Load Ausdiab lipid names
Ausdiab_lipid_names <- read.table("Ausdiab_lipid_names.tsv", sep="\t")$V1


##############################################################################
##
##  Create LIPID database
##
##############################################################################

# Create a data frame with completely simulated data
data_LIPID <- data.frame(
    ID = seq(1, n),
    pid = seq(1, n),
    period = sample(c(0, 12), n, replace = TRUE),
    TREAT = sample(c(0, 1), n, replace = TRUE),
    age = sample(20:80, n, replace = TRUE),
    sexn = sample(c("Male", "Female"), n, replace = TRUE),
    bmi = rnorm(n, mean = 25, sd = 5)
)

# Load LIPID lipid names
LIPID_lipid_names <- read.table("LIPID_lipid_names.tsv", sep="\t")$V1


##############################################################################
##
##  Load Supplementary table 2
##
##############################################################################

# Create a full list of lipids to simulate
lipids_to_simulate <- unique(c(Ausdiab_lipid_names, LIPID_lipid_names))

# Number of variables (lipid species)
p <- length(lipids_to_simulate)

## Get the size of blocks
n_blocks <- ceiling(p / 50)

# Create a block diagonal matrix with high correlations in blocks
blocks <- lapply(rep(n_blocks, 50), function(size) {
    
    # Simulate the lower triangle
    lower_tri <- runif(size * (size - 1) / 2, min = 0.4, max = 0.8)
    
    # Generate a symmetric matrix
    symmetric_matrix <- diag(0.5, size, size)
    
    # Fill in the lower triangle
    symmetric_matrix[lower.tri(symmetric_matrix)] <- lower_tri
    
    # Fill in the upper triangle
    symmetric_matrix <- symmetric_matrix + t(symmetric_matrix)
    
    # Return matrix
    return(symmetric_matrix)
})

# Fill out the correlation matrix
corr_matrix <- do.call(Matrix::bdiag, blocks)
corr_matrix <- as.matrix(corr_matrix[seq(p), seq(p)])

# Ensure the matrix is positive definite (required for mvrnorm)
eigenvalues <- eigen(corr_matrix, symmetric = TRUE)$values
min_eigenvalue <- min(eigenvalues)
if(min_eigenvalue <= 0) {
    corr_matrix <- corr_matrix - (min_eigenvalue - 0.01) * diag(p)
    corr_matrix <- cov2cor(corr_matrix)
}

# Mean vector (assuming mean of 0 for simplicity)
mu <- rep(0, p)

# Generate data from multivariate normal distribution with specified correlation matrix
simulated_data <- mvrnorm(n * 2, mu = mu, Sigma = corr_matrix)

# Add lipid names
colnames(simulated_data) <- lipids_to_simulate

# Split into two data sets
simulated_data_Ausdiab <- simulated_data[1:n, Ausdiab_lipid_names]
simulated_data_LIPID <- simulated_data[(n + 1):(2 * n), LIPID_lipid_names]

# Exponentiate the data
simulated_data_Ausdiab <- exp(simulated_data_Ausdiab)
simulated_data_LIPID <- exp(simulated_data_LIPID)

# Add lipid data to simulated clinical data
simulated_data_Ausdiab <- cbind("Not real data", data_Ausdiab, simulated_data_Ausdiab)
simulated_data_LIPID <- cbind("Not real data", data_LIPID, simulated_data_LIPID)


##############################################################################
##
##  Save the Simulated data
##
##############################################################################

# Save the Simulated Ausdiab data to a CSV file (replaces 2021_12_16 - Ausdiab lipidomics and clinical data combined.csv in scripts 1-6)
write.csv(simulated_data_Ausdiab, "Simulated_Dataset_1.csv", row.names = FALSE)

# Save the Simulated LIPID data to a CSV file (replaces 2017_09_12 - Lipid Trial Database - Mapped to Ausdiab.csv in scripts 1-6)
write.csv(simulated_data_LIPID, "Simulated_Dataset_2.csv", row.names = FALSE)

