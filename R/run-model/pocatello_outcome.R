#!/usr/bin/env Rscript

max_cores <- 18

# INTERVENTION IMPACT (cluster_outcome)
# Estimate intervention impact with posterior samples
# Run posterior LHS samples through the model
# Include intervention scenarios of interest 
# test pop:   1 = all agyw, 2 = sa agyw, 3 = pregnant, 4 = fsw, 5 = male non targeted, 6 = male high risk
# test strat: h = health service screening, y = syndromic management recipients

run_partition_outcome <- function(indices_to_run){
  
  source("cluster_outcome.R")
  library(tidyverse)
  library(odin)
  library(pkgload)
  library(data.table)
  
  launch_outcome(df_param_base, df_excel_init, df_ent, df_preg, 
                 post_samples, 
                 run_model, outcome_diag, outcome_screen, 
                 prep_param, prep_init, tt,
                 test_pop = c(2,3,4,5,6), test_num = c(0,10), test_strat = c("h","y"),
                 indices_to_run = indices_to_run)
  
}

# Split rows across cores
post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)
parallel_indices <- parallel::splitIndices(nrow(post_samples), max_cores)

cl <- parallel::makeCluster(max_cores)

id_result <- parallel::clusterApply(cl, parallel_indices, run_partition_outcome)

# Merge result_list from each cluster in the correct order
df_merge_result <- unlist(lapply(id_result, function(result) result$result_list), recursive = FALSE)

df_bind_result <- dplyr::bind_rows(df_merge_result) |>
  dplyr::mutate(test_pop = factor(test_pop, levels = c(0,1,2,3,4,5,6),
                                  labels = c("none","agyw_all","agyw_sa","preg","fsw","male_all","male_h")))

#write.csv(df_bind_result, file = "./posterior/post_intervention_constrained.csv", row.names = FALSE)
write.csv(df_bind_result, file = "./posterior/post_intervention_unrestricted.csv", row.names = FALSE)

parallel::stopCluster(cl)