#!/usr/bin/env Rscript

v <- "v1"

max_cores <- 18

# PREV ESTIMATES (cluster_launch)
# Estimate prevalence with prior samples
# i.e. run full data frame of prior LHS samples through model to calculate prevalence, number sex partners, kappa_p, and condom use
run_partition <- function(indices_to_run){
  
  source("cluster_launch.R")
  library(tidyverse)
  library(odin)
  library(pkgload)
  
  launch_process(df_param_base, df_excel_init, df_ent, df_preg, 
                 prior_samples, run_model, prep_param, prep_init,
                 tt,
                 indices_to_run = indices_to_run)
  
}

# Split rows across cores
prior_samples <- read.csv("./prior/prior_samples.csv"), check.names = FALSE)
parallel_indices <- parallel::splitIndices(nrow(lhs_samples), max_cores)

cl <- parallel::makeCluster(max_cores)

id_result <- parallel::clusterApply(cl, parallel_indices, run_partition)

# Mege df_prev from across cores
merge_prev <- function(result, index) {
  node_indices <- parallel_indices[[index]]
  result$df_prev[node_indices]
}

df_merge_prev <- unlist(mapply(merge_prev, id_result, seq_along(id_result)), recursive = FALSE)
df_bind_prev <- dplyr::bind_rows(df_merge_prev)

write.csv(df_bind_prev, file = "./prior/prior_prev.csv" row.names = FALSE)

# MERGE df_c 
merge_c <- function(result, index) {
  node_indices <- parallel_indices[[index]]
  result$df_c[node_indices]
}

df_merge_c <- unlist(mapply(merge_c, id_result, seq_along(id_result)), recursive = FALSE)
df_bind_c <- dplyr::bind_rows(df_merge_c)

write.csv(df_bind_c, file = "./prior/prior_c.csv", row.names = FALSE)


# MERGE df_kappa 
merge_kappa <- function(result, index) {
  node_indices <- parallel_indices[[index]]
  result$df_kappa[node_indices]
}

df_merge_kappa <- unlist(mapply(merge_kappa, id_result, seq_along(id_result)), recursive = FALSE)
df_bind_kappa <- dplyr::bind_rows(df_merge_kappa)

write.csv(df_bind_kappa, file = "./prior/prior_kappa.csv", row.names = FALSE)