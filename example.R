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

## Run Bisect to estimate cell type proportions for the individuals in the dataset.
P <- run_EM(R = R, Y = Y, Pi = Pi, estimate_Pi = F)


## Compare to real cell type composition --------------------------------------------------------------
baseline <- read_csv('data/baseline_proportions.csv')

visualization_result <- get_visualization_dataframe(P, baseline)

## plot a scatter plot of true cell types vs estimated.  Looks pretty good!
visualization_result %>% ggplot(aes(truth, estimate, color=cell_type)) + geom_point(size=0.2, alpha = 0.4) + 
  geom_abline(intercept = 0, slope = 1) + xlab("True Cell Proportion") + ylab("Estimated Cell Proportion") + 
  guides(colour = guide_legend(override.aes = list(size=10))) + scale_color_discrete(name = "Cell Type")
