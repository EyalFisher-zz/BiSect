library(tidyverse)
source('EM.R')
source('visualization.R')

## Load the methylation reference, that contains methylation probabilities per cell type in 241 chosen sites.
reference <- read_csv('data/methylation_reference.csv')

## Load the matrices of methylation and total reads from the same sites.
##Notice that the sites in all 3 matrices (reference, methylation, reads), need to appear in the same order(!).
methylation <- read_csv('data/GSE402079_methylation_matrix.csv')
reads <- read_csv('data/GSE402079_reads_matrix.csv')

## The matrices need to be in the form of a matrix of individuals (rows) X sites (sites).
Y <- t(as.matrix(methylation))
R <- t(as.matrix(reads))

## The column of the IDs needs to be removed from Pi, as the matrix will be used for computations
Pi <- as.matrix(reference[,-1])

## the hyper parameters for the dirichlet distribution. Was estimated by fitting a Dirichlet distribution to cell counts data.
alpha=c(2.5392, 1.7934, 0.7240, 0.7404, 1.8439, 15.0727)

initial_P <- initialize_P(nrow(Y), ncol(Pi), alpha)

## Run Bisect to estimate cell type proportions for the individuals in the dataset.
P <- run_EM(R = R, Y = Y, Pi = Pi, P = initial_P ,estimate_Pi = F, alpha = alpha, iterations = 200)

## Compare to real cell type composition --------------------------------------------------------------
baseline <- read_csv('data/baseline_proportions.csv')

visualization_result <- get_visualization_dataframe(P, baseline)

# plot a scatter plot of true cell types vs estimated.  Looks pretty good!
visualization_result %>% ggplot(aes(truth, estimate, color=cell_type)) + geom_point(size=0.2, alpha = 0.4) + 
  geom_abline(intercept = 0, slope = 1) + xlab("True Cell Proportion") + ylab("Estimated Cell Proportion") + 
  guides(colour = guide_legend(override.aes = list(size=10))) + scale_color_discrete(name = "Cell Type")


## Use semi-supervised mode ---------------------------------------------------------------------------
set.seed(1234)
# Choose 50 random individuals with known cell type composition
n_known_samples <- 50
known_samples_indices <- sample.int(nrow(baseline), size = n_known_samples)   
known_samples <- as.matrix(baseline[known_samples_indices, ])

# Fit a dirichlet distirbutio nto the known samples to use as a prior
fit_dirichlet <- sirt::dirichlet.mle(as.matrix(known_samples))
alpha <- fit_dirichlet$alpha

# Organize all the matrices such that the known samples are the first 50 rows.
Y_known <- Y[known_samples_indices, ]
Y_unknown <-Y[-known_samples_indices, ]
R_known <- R[known_samples_indices, ]
R_unknown <- R[-known_samples_indices, ]

all_Y <- rbind(Y_known, Y_unknown)
all_R <- rbind(R_known, R_unknown)

# generate an initial P matrix for the unknown individuals
initial_P <- initialize_P(n_individuals = nrow(R_unknown), n_cell_types = length(alpha), alpha)

# attach the known proportions so that they are the first 50 rows.
P <- rbind(known_samples, initial_P)

# initialize unknown Pi (change for methylation in each cell type, in each site) to random values
Pi <- initialize_Pi(n_sites = ncol(R), n_cell_types = length(alpha))

# Run Bisect, making sure to supply the number of known individuals.
P <- run_EM(R = all_R, Y = all_Y, Pi = Pi, P = P, n_known_samples = n_known_samples, estimate_Pi = T, alpha = alpha, iterations = 200)

## Compare to real cell composition
baseline <- read_csv('data/baseline_proportions.csv')

baseline_unknown <- baseline[-known_samples_indices, ]
all_baseline <- rbind(known_samples, baseline_unknown)

visualization_result <- get_visualization_dataframe(P, all_baseline)

## plot a scatter plot of true cell types vs estimated.  Looks pretty good!
visualization_result %>% ggplot(aes(truth, estimate, color=cell_type)) + geom_point(size=0.2, alpha = 0.4) + 
  geom_abline(intercept = 0, slope = 1) + xlab("True Cell Proportion") + ylab("Estimated Cell Proportion") + 
  guides(colour = guide_legend(override.aes = list(size=10))) + scale_color_discrete(name = "Cell Type")
