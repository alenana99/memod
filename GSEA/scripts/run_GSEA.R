#!/usr/bin/env Rscript

# Definition of run_GSEA function
run_GSEA <- function(list_S, L) {
  # Esegui ES_significance2 su ciascun data frame
  results <- lapply(list_S, function(data_list) {
    lapply(data_list, function(df) {
      ES_significance2(df, L)
    })
  })
  
  # Extract NES values from each data frame
  NES_values <- lapply(results, function(sublist) {
    sapply(sublist, function(x) x$NES)
  })
  
  # Extract NES_pi values from each data frame
  NES_pi_values <- lapply(results, function(sublist) {
    do.call(c, lapply(sublist, function(x) x$NES_pi))
  })
  
  # Calculate FDR for each sublist
  FDR_results <- lapply(results, function(sublist) {
    NES_sub <- sapply(sublist, function(x) x$NES)
    NES_pi_sub <- do.call(c, lapply(sublist, function(x) x$NES_pi))
    FDR_wrapper(NES_sub, NES_pi_sub)
  })
  
  return(list(NES_values = NES_values, NES_pi_values = NES_pi_values, FDR_results = FDR_results))
}

args <- commandArgs(trailingOnly = TRUE)
list_S_path <- args[1]
L_path <- args[2]

list_S <- readRDS(list_S_path)
L <- read.csv(L_path, header = TRUE)

# Esegui GSEA
gsea_results <- run_GSEA(list_S, L)
