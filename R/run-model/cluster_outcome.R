# INTERVENTION IMPACT (cluster_outcome)
# Estimate intervention impact with posterior samples
# Run posterior LHS samples through the model
# Include intervention scenarios of interest 
# test pop:   1 = all agyw, 2 = sa agyw, 3 = pregnant, 4 = fsw, 5 = male non targeted, 6 = male high risk
# test strat: h = health service screening, y = syndromic management recipients

source("prep_input.R") # for prep_param() and prep_init()
df_excel_param <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "param")
df_excel_init <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "init")
df_ent <- read.csv("./parameters/entrants.csv")
df_preg <- read.csv("./parameters/pregnancy.csv")
df_param_base <- prep_param(df_excel_param, "baseline") 
post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)
#post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)[c(1),]
tt <- c(0,10,20,40,50,54,55,56,57,58,59,60)
#tt <- c(0:60)

# Run odin model
run_model <- function(df_input, tt) {
  
  #gen <- odin::odin("model_odin_constrained.R")
  gen <- odin::odin("model_odin_unrestricted.R")
  
  mod <- gen$new(user = df_input)
  df <- mod$run(tt) |>
    as.data.frame()  |>
    pivot_longer(!t) |> 
    mutate(var = sub("\\[.*", "", name)) |>
    filter(var %in% c("N", "NS", "I",
                      "X", "Y", "Z",
                      "My", "Ty", "Ry",
                      "Th_symp", "Th_asymp",
                      "Rh_symp", "Rh_asymp",
                      "N_hc", 
                      #"tau_hc",
                      "n_I", "n_X", "n_Y", "n_Z",
                      "n_Ty", "n_My", "n_Ry",
                      "n_Th_symp", "n_Th_asymp",
                      "n_Rh_symp", "n_Rh_asymp",
                      "n_symp", 
                      #"delta",
                      "n_test_h", "n_test_y", "n_test_o",
                      "n_elig_h", "n_elig_y", "n_elig_o",
                      "py_X", "py_Y", "py_Z",
                      "py_My", "py_Ty", "py_Ry",
                      "py_Th_symp", "py_Th_asymp",
                      "py_Rh_symp", "py_Rh_asymp")) |>
    mutate(
      strat = stringr::str_extract(name, "\\[.*\\]"),
      sex = if_else(nchar(strat) > 0, substring(strat, 2, 2), NA_character_),
      risk = if_else(nchar(strat) > 3, substring(strat, 4, 4), NA_character_),
      age = if_else(nchar(strat) > 5, substring(strat, 6, 6), NA_character_),
      gest = if_else(nchar(strat) > 7, substring(strat, 8, 8), NA_character_),
      sex2 = if_else(nchar(strat) > 9, substring(strat, 10, 10), NA_character_),
      risk2 = if_else(nchar(strat) > 9, substring(strat, 12, 12), NA_character_),
      age2 = if_else(nchar(strat) > 9, substring(strat, 14, 14), NA_character_),
      gest2 = if_else(nchar(strat) > 9, substring(strat, 16, 16), NA_character_),
      year = t + 1970) |>
    select(t, year, var, sex:gest2, value)
  
}

# Calculate outcomes
outcome_calc <- function(df_result) {

  # Convert to data.table for faster operations
  df_result <- as.data.table(df_result)
  
  # Define columns to sum
  sum_cols <- c("N", "NS", "I", "X", "Y", "Z", 
                "Th_symp", "Th_asymp", "Rh_symp", "Rh_asymp",
                "Ty", "My", "Ry", 
                "n_I", "n_X", "n_Y", "n_Z", 
                "N_hc", 
                "n_Th_symp", "n_Th_asymp", "n_Rh_symp", "n_Rh_asymp",
                "n_Ty", "n_My", "n_Ry", 
                "n_symp",
                "n_elig_h", "n_elig_y", "n_elig_o", 
                "n_test_h", "n_test_y", "n_test_o",
                "py_X", "py_Y", "py_Z", 
                "py_Th_symp", "py_Th_asymp", "py_Rh_symp", "py_Rh_asymp",
                "py_My", "py_Ty", "py_Ry")
  
  # Filter and summarize in one operation
  df_result[year %in% 2024:2030, 
            lapply(.SD, sum), 
            by = .(lhs_num, year, test_pop, test_strat, test_num, sex),
            .SDcols = sum_cols]
}

# Outcome loop for viable parameter samples
launch_outcome <- function(df_param_base, df_excel_init, df_ent, df_preg, 
                           post_samples, 
                           run_model, outcome_diag, outcome_screen, 
                           prep_param, prep_init, tt,
                           # test strategies
                           test_pop, test_num, test_strat,
                           indices_to_run) {
  
  # Pre-calculate total number of scenarios
  n_scenarios <- length(test_pop) * length(test_num) * length(test_strat)
  total_scenarios <- length(indices_to_run) * n_scenarios
  
  # Pre-allocate outcome list
  #outcome_list <- vector("list", total_scenarios)
  result_list <- vector("list", total_scenarios)
  
  # Create all possible scenario combinations upfront
  scenarios <- expand.grid(
    test_strat = test_strat,
    test_pop = test_pop,
    test_num = test_num,
    stringsAsFactors = FALSE
  )
  
  # Counter for outcome_list
  counter <- 0
  
  # Loop through each row
  for (i in indices_to_run) {
    
    # fresh parameter set
    df_param_new <- df_param_base 
    
    # drop second value for pi and sigma
    df_param_new$pi <- df_param_base$pi[1]
    df_param_new$sigma <- df_param_base$sigma[1]
    
    # Prepare parameter set for row of values
    for (j in 1:(ncol(post_samples)-1)) {    # ignore lh_num at end
      
      # set the variable, stratum and value to update in df_param_new
      var <- sub("-.*", "", colnames(post_samples)[[j]]) 
      strat <- as.numeric(sub(".*-", "", colnames(post_samples)[[j]]))
      value <- post_samples[i,j]
      
      # update the variable
      df_param_new[[var]][[strat]] <- value
    }
    
    # Update array calculations with new parameter set
    df_param_update <- prep_param(df_param_new, "update")
    
    # Re-calculate initial conditions
    df_init <- prep_init(df_excel_init, df_param_update)
    
    # Prepare full input parameter set
    df_input <- c(df_param_update, df_init, df_ent, df_preg,
                  list(test_t = c(0,55,60), test_h = c(0,0,0), test_y = c(0,0,0), test_pop = 0))
    
    # Loop through scenarios
    for (s in 1:nrow(scenarios)) {
      
      # Create clean parameter list
      df_clean <- df_input
      
      # Update test parameters for this scenario
      df_clean$test_pop <- scenarios$test_pop[s]
      strat_var <- paste0("test_", scenarios$test_strat[s])
      df_clean[[strat_var]][2:3] <- scenarios$test_num[s]
      
      # Run model
      result <- run_model(df_clean, tt) |>
        pivot_wider(names_from = var, values_from = value) |>
        mutate(lhs_num = post_samples[i,"lhs_num"], 
               test_pop = scenarios$test_pop[s], 
               test_strat = scenarios$test_strat[s], 
               test_num = scenarios$test_num[s])
      
      # Calculate outcomes
      counter <- counter + 1
      #outcome_list[[counter]] <- outcome_calc(result)
      result_list[[counter]] <- result
      
    }
  }
  
  return(list(result_list = result_list))
  #return(list(outcome_list = outcome_list))
}
