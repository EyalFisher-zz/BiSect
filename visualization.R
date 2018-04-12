
get_visualization_dataframe <- function(bisect_results, true_cell_counts)
{
  estimates_bin <- as.tibble(bisect_results)
  true_cell_counts <- as.tibble(true_cell_counts)
  
  colnames(estimates_bin) <- c('CD4', 'CD8', 'mono', 'Bcells', 'NK', 'gran')
  colnames(true_cell_counts) <- c('CD4', 'CD8', 'mono', 'Bcells', 'NK', 'gran')
  
  gathered_estimates_bin <- estimates_bin %>% gather('CD4', 'CD8', 'mono', 'Bcells', 'NK', 'gran', key='cell_type', value='estimate_norm')
  gathered_truth <- true_cell_counts %>% gather('CD4', 'CD8', 'mono', 'Bcells', 'NK', 'gran', key='cell_type', value='truth')
  
  gathered_estimates_bin <- gathered_estimates_bin %>% mutate(method='bin')
  colnames(gathered_estimates_bin) <- c('cell_type', 'estimate', 'method')
  
  estimates <- rbind(gathered_estimates_bin)
  truth <- rbind(gathered_truth, gathered_truth)
  
  results <- cbind(truth, select(estimates, 'estimate', 'method'))
  
  return(results)
}