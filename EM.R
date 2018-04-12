E_1_step <- function(P, Pi)
{
  N <- nrow(P)
  J <- nrow(Pi)
  K <- ncol(P)
  
  E <- array(dim=c(N,J,K))
  
  constants <- P %*% t(Pi)
  
  results <- c()
  for(k in 1:K)
  {
    E[,,k] <- outer(P[,k], Pi[,k])/constants
  }
  
  return(E)
}


E_0_step <- function(P, Pi)
{
  N <- nrow(P)
  J <- nrow(Pi)
  K <- ncol(P)
  
  E <- array(dim=c(N,J,K))
  
  constants <- P %*% (1-t(Pi))
  
  results <- c()
  for(k in 1:K)
  {
    E[,,k] <- outer(P[,k], 1-Pi[,k])/constants
  }
  
  return(E)
}

Maximize_P <- function(E_0, E_1, R, Y, reads_i, alpha)
{
  K <- dim(E_0)[3]
  constants <- reads_i + sum(alpha) - K
  
  N <- nrow(R)
  
  result <- matrix(nrow=N, ncol=K)
  
  for(k in 1:K)
  {
    nominator <- rowSums(E_1[,,k] * Y) + rowSums(E_0[,,k] * (R-Y)) + alpha[k] - 1
    result[,k] <- nominator/constants
  }
  
  return(result)
}

Maximize_Pi <- function(E_0, E_1, Y, R)
{
  J <- ncol(Y)
  K <- dim(E_0)[3]
  result <- matrix(nrow=J, ncol=K)
  for(k in 1:K)
  {
    nominator <- colSums(Y * E_1[,,k])
    denominator <- colSums(Y * E_1[,,k]) + colSums((R- Y) * E_0[,,k])
    
    result[,k] <- nominator/denominator
  }
  
  return(result)
}


## R - a matrix of individuals (rows) on sites (columns), containing the total number of reads for each site, in each individual.
## Y - a matrix of individuals (rows) on sites (columns), containing the number of methylated reads for each site, in each individual.
## Pi - A matrix of sites (rows) on cell types (columns), containing the probability for methylation in each site, in each cell type. 
## If not known, initialize to a random value between 0 and 1.
## alpha - a list of length K (amount of cell types), containing the hyper-paramters of the dirichlet prior. If None, a non-informative prior is used.
run_EM <- function(R, Y, Pi, iterations=200, estimate_Pi=F ,alpha = NA)
{
  n <- nrow(R)
  K <- ncol(Pi)
  J <- ncol(R)
  
  minimum_cell_proportion <- 0.001
  if(length(alpha) == 1)
  {
    if(is.na(alpha))
    {
      alpha <- rep(1, K) 
    }
  }
  
  ## Initialize P to the mean cell composition in the population (according to prior)
  P <- rep(alpha/sum(alpha), n)
  P <- matrix(P, nrow = n)
  
  reads_i <- rowSums(R)
  
  for(i in 1:iterations)
  {
    E_0 <- E_0_step(P, Pi)
    E_1 <- E_1_step(P, Pi)
    
    P <- Maximize_P(E_0, E_1, R, Y, reads_i, alpha=alpha)
    
    ## make sure that P doesn't leave allowed range to prevent numeric errors.
    P[P < 0] <- minimum_cell_proportion 
    P[P > 1] <- 1-minimum_cell_proportion
    
    
    if(estimate_Pi) {
      Pi <- Maximize_Pi(E_0, E_1, Y, R)
    }
    
    print(i)
    
  }
  
  return(P)
}