library(lhs)
library(tidyverse)

# PREP PRIOR LHS SAMPLES ----

# < Sample set up ----

# List of variables to sample
df_lhs <- readxl::read_excel("parameters.xlsx", sheet = "param") |>
  select(variable, var, var_strat, min, max) |>
  filter(!is.na(var_strat)) 

# df_lhs <- df_lhs |> 
#   filter(var %in% c("kappa","phi","pi",
#                     "c_low","c_high_ratio",
#                     "chi_low", "chi_high_ratio",
#                     #"e",
#                     #"upsilon",
#                     #"theta_c", "theta_chi",
#                     #"n_l",
#                     "gamma_low", "gamma_sex_ratio", "gamma_high_ratio",
#                     "tau_m_low", "tau_m_sex_ratio", #"tau_m_high_ratio",
#                     #"ap",
#                     "zeta"
#                     ))

# Number of samples per variable
# Set number of simulations (or runs=N) to 4/3 times greater than number of 
# uncertain parameters k, i.e. N > 4/3k [Blower and Dowlatabadi 1994].
n_samples <- 50000

# Create matrix for the samples
# column corresponds to each row in df_lhs
lhs_samples <- randomLHS(n_samples, length(df_lhs$max))

# Prepare to scale samples to parameter ranges
scaled_samples <- as.data.frame(lhs_samples)
colnames(scaled_samples) <- df_lhs$var_strat

# Extract min and max values for scaling
min_vals <- df_lhs$min
max_vals <- df_lhs$max

# Scale samples
# Function to scale samples between min and max values
scale_samples <- function(samples, min_vals, max_vals) {
  scaled <- t(apply(samples, 1, function(x) min_vals + x * (max_vals - min_vals)))
  return(as.data.frame(scaled))
}

final_samples <- scale_samples(scaled_samples, min_vals, max_vals)
# Note that don't need expand_grid() on final_samples
# LHS already accounts for having multiple variables. It ensures these variables
# are evenly distributed among each other.

# < Check for errors in parameter definitions ----

# Check whether any variable and relative rate combinations cause errors
# i.e. probability > 1 or probability < 0
# drop those

source("prep_input.R")
df_excel_param <-  readxl::read_excel("./parameters/parameters.xlsx", sheet = "param")
df_param_base <- prep_param(df_excel_param, "baseline")

# Loop through to update parameters for each row of the lhs samples

param_check <- list()

for (i in 1:nrow(final_samples)){
  
  df_param_new <- df_param_base # fresh parameter set
    
    for (j in 1:ncol(final_samples)) {    
      
      # set the variable, stratum and value to update in df_param_new
      var <- sub("-.*", "", colnames(final_samples)[[j]]) 
      strat <- as.numeric(sub(".*-", "", colnames(final_samples)[[j]]))
      value <- final_samples[i,j]
      
      # update the variable
      df_param_new[[var]][[strat]] <- value
    }
    
  # Update array calculations with new parameter set
  df_param_update <- prep_param(df_param_new, "update")

  param_check[[i]] <- df_param_update
    
}

# Check which variables cause errors
error_check <- lapply(param_check, function(x) x$error_check)

error_list <- cbind(error_check) |> as.data.frame() |> mutate(lhs_num = row_number()) |>
  filter(!is.na(error_check))

error_list |> count(error_check)

# Review gamma errors
gamma_list <- cbind(error_check) |> as.data.frame() |> mutate(lhs_num = row_number()) |>
  filter(error_check %in% c("gamma","chi,gamma")) |>
  pull(lhs_num)

lapply(param_check, function(x) x$gamma)[gamma_list]

# Review chi errors
chi_list <- cbind(error_check) |> as.data.frame() |> mutate(lhs_num = row_number()) |>
  filter(error_check %in% c("chi","chi,gamma")) |>
  pull(lhs_num)

lapply(param_check, function(x) x$chi)[chi_list]

# gamma_check <- lapply(param_check, function(x) x$gamma)
# chi_check <- lapply(param_check, function(x) x$chi)
# chi_check_fsw <- lapply(param_check, function(x) x$chi[2,3,1,1])
# c_check <- lapply(param_check, function(x) x$c)
# c_check_fsw <- lapply(param_check, function(x) x$c[2,3,1,1])

# Drop rows with errors from final_samples (as long as not too many)
final_samples <- final_samples |>
  filter(!row_number() %in% (error_list |> pull(lhs_num)))

# I get 14.9% of rows that have errors
# Therefore need initial sample as:
# 10000 = x - x*0.149
# 10000 = x*(1 - 0.149)
# x = 10000/(1-0.149)

write.csv(final_samples, file = "./prior/prior_samples.csv", row.names = FALSE)