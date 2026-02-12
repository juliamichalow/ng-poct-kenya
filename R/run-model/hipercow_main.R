# SET UP HIPERCOW ----
library(hipercow)
#hipercow_init(driver = "windows")
windows_check()

hipercow_provision() # install packages from pkgdepends.txt
# hipercow_provision("pkgdepends", refs = "cran::odin") # manually install
# hipercow_provision("pkgdepends", refs = "cran::tidyverse") # manually install

hipercow_provision_list()
hipercow_provision_check()
hipercow_provision_compare()

max_cores <- 32

resources <- hipercow_resources(cores = max_cores)
options(hipercow.max_size_local = Inf)

hipercow_environment_create(packages = c("odin","tidyverse","pkgload"),
                            sources = c("cluster_launch.R",
                                        "prep_input-3.R"))

# hipercow_environment_create(packages = c("odin","tidyverse","pkgload"),
#                            sources = c("cluster_outcome.R",
#                                        "prep_input-3.R"))
#                            #globals = TRUE)

hipercow_configuration()

# PREV ESTIMATES (cluster_launch) ----
# Estimate prevalence with prior samples
# i.e. run full data frame with LHS variable samples through the model to calculate prevalence 
## < Run in parallel ----
run_partition <- function(indices_to_run){
  
  # partition the function based on indices_to_run
  # i.e. divide number of rows in lhs_samples across number of cores
  # run the loop in parallel across cores
  
  launch_process(df_param_base, df_excel_init, df_ent, df_preg, 
         lhs_samples, run_model, prep_param, prep_init,
         tt,
         indices_to_run = indices_to_run)
         
}

# Split rows across cores
lhs_samples <- read.csv("lhs_samples.csv", check.names = FALSE)
parallel_indices <- parallel::splitIndices(nrow(lhs_samples), max_cores)

id <- task_create_expr(
  parallel::clusterApply(NULL, parallel_indices, run_partition),
  parallel = hipercow_parallel("parallel"),
  resources = resources)

write.csv(id, file = "id_process.csv")
# id <- read.csv("id_process.csv")[[2]]

task_status(id)
task_log_show(id)

id_result <- task_result(id)

## < Process results ----
# Merge results that are split across cores into single list

library(tidyverse)

# MERGE df_prev  
merge_prev <- function(result, index) {
  node_indices <- parallel_indices[[index]]
  result$df_prev[node_indices]
}

df_merge_prev <- unlist(mapply(merge_prev, id_result, seq_along(id_result)), recursive = FALSE)
df_bind_prev <- bind_rows(df_merge_prev)

write.csv(df_bind_prev, 
          file = "C:/Users/jhm21/OneDrive - Imperial College London/Analysis/Objective 4/Model - srag/ng-transmission-srag-v2/lhs-3/lhs_10k/v4/lhs_prevalence.csv", 
          row.names = FALSE)


# MERGE df_c 
merge_c <- function(result, index) {
  node_indices <- parallel_indices[[index]]
  result$df_c[node_indices]
}

df_merge_c <- unlist(mapply(merge_c, id_result, seq_along(id_result)), recursive = FALSE)
df_bind_c <- bind_rows(df_merge_c)

write.csv(df_bind_c, 
          file = "C:/Users/jhm21/OneDrive - Imperial College London/Analysis/Objective 4/Model - srag/ng-transmission-srag-v2/lhs-3/lhs_10k/v4/lhs_c_ad.csv", 
          row.names = FALSE)


# MERGE df_kappa 
merge_kappa <- function(result, index) {
  node_indices <- parallel_indices[[index]]
  result$df_kappa[node_indices]
}

df_merge_kappa <- unlist(mapply(merge_kappa, id_result, seq_along(id_result)), recursive = FALSE)
df_bind_kappa <- bind_rows(df_merge_kappa)

write.csv(df_bind_kappa, 
          file = "C:/Users/jhm21/OneDrive - Imperial College London/Analysis/Objective 4/Model - srag/ng-transmission-srag-v2/lhs-3/lhs_10k/v4/lhs_kappa.csv", 
          row.names = FALSE)


# MEREG df_result (odin model results)
merge_result <- function(result, index) {
  node_indices <- parallel_indices[[index]]
  result$df_result[node_indices]
}

df_merge_result <- unlist(mapply(merge_result, id_result, seq_along(id_result)), recursive = FALSE)

# Save result for each model run
library(data.table)
for (i in seq_along(df_merge_result)) {
  fwrite(df_merge_result[[i]], 
         paste0("C:/Users/jhm21/OneDrive - Imperial College London/Analysis/Objective 4/Model - srag/ng-transmission-srag-v2/lhs-3/lhs_10k/result_per_run/result_",i,".csv"), 
         append = TRUE)
}



# INTEVRENTION IMPACT (cluster_outcome) ----
# Estimate intervention impact with posterior samples
# Run posterior dataframe with viable LHS variable samples through the model
# Include intervention scenarios of interest 
# test pop:   1 = all agyw, 2 = sa agyw, 3 = pregnant, 4 = fsw, 5 = male non targeted, 6 = male high risk
# test strat: h = health service screening, y = syndromic management recipients

## < Run in parallel ----
run_partition_outcome <- function(indices_to_run){
  
  # partition the function based on indices_to_run
  # i.e. divide number of rows in lhs_samples across number of cores
  # run the loop in parallel across cores
  
  launch_outcome(df_param_base, df_excel_init, df_ent, df_preg, 
                 lhs_samples_post |> select(-num), 
                 run_model, prep_param, prep_init, tt,
                 test_pop = c(2,3,4,5,6), test_num = c(0,150000,450000), test_strat = c("h","y"),
                 indices_to_run = indices_to_run)
  
}

# Split rows across cores
lhs_samples_post <- read.csv("lhs_samples_post.csv", check.names = FALSE)
parallel_indices_outcome <- parallel::splitIndices(nrow(lhs_samples_post), max_cores)

id <- task_create_expr(
  parallel::clusterApply(NULL, parallel_indices_outcome, run_partition_outcome),
  parallel = hipercow_parallel("parallel"),
  resources = resources)

# id is for 10 k v1
write.csv(id, file = "id_outcome.csv")
#id <- read.csv("id_outcome.csv")[[2]]
task_status(id)
task_log_show(id)
id_result <- task_result(id)

## < Process results ----

# merge result_list from each cluster in the correct order
df_merge_result <- unlist(lapply(id_result, function(result) result$result_list), recursive = FALSE)

# merge outcome_list from each cluster in the correct order
df_merge_outcome <- unlist(lapply(id_result, function(result) result$outcome_list), recursive = FALSE)

library(tidyverse)
df_bind_outcome <- bind_rows(df_merge_outcome) |>
  mutate(test_pop = factor(test_pop, levels = c(0,1,2,3,4,5,6),
                           labels = c("none","agyw_all","agyw_sa","preg","fsw","male_all","male_h"))) 

write.csv(df_bind_outcome, 
          file = "C:/Users/jhm21/OneDrive - Imperial College London/Analysis/Objective 4/Model - srag/ng-transmission-srag-v2/lhs-3/lhs_10k/v1/lhs_outcome.csv", 
          row.names = FALSE)
