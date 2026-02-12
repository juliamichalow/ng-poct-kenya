# PREV ESTIMATES (cluster_launch)
# Estimate prevalence with prior samples
# i.e. run full data frame of prior LHS samples through model to calculate prevalence, number sex partners, kappa_p, and condom use

source("prep_input.R") # for prep_param() and prep_init()
df_excel_param <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "param")
df_excel_init <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "init")
df_ent <- read.csv("./parameters/entrants.csv")
df_preg <- read.csv("./parameters/pregnancy.csv")
df_param_base <- prep_param(df_excel_param, "baseline")
prior_samples <- read.csv("./prior/prior_samples.csv", check.names = FALSE)
tt <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
# tt <- c(0, 10, 20, 30, 40, 50, 60)

run_model <- function(df_input, tt) {
  
  gen <- odin::odin("model_odin_constrained.R")
  
  mod <- gen$new(user = df_input)
  df <- mod$run(tt) |>
    as.data.frame()  |>
    pivot_longer(!t) |> 
    mutate(
      var = sub("\\[.*", "", name),
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
  
  return(df)
  
}

launch_process <- function(df_param_base, df_excel_init, df_ent, df_preg, 
                   prior_samples, run_model, prep_param, prep_init,
                   tt, indices_to_run) {
  
  # Empty list for each model's result and calculated prevalence
  # df_result <- list()
  df_prev <- list()
  df_c <- list()
  df_kappa <- list()
  
  # Loop through each row
  for (i in indices_to_run) {
    
    df_param_new <- df_param_base # fresh parameter set
    # drop second value for pi and sigma
    df_param_new$pi <- df_param_base$pi[1]
    df_param_new$sigma <- df_param_base$sigma[1]
    
    for (j in 1:ncol(prior_samples)) {    
      
      # set the variable, stratum and value to update in df_param_new
      var <- sub("-.*", "", colnames(prior_samples)[[j]]) 
      strat <- as.numeric(sub(".*-", "", colnames(prior_samples)[[j]]))
      value <- prior_samples[i,j]
      
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
    
    # Run model
    result <- run_model(df_input, tt)
    
    # Store model result
    # df_result[[i]] <- result
    
    # Calculate prevalence
    df_prev[[i]] <- result |>
     filter(var %in% c("I", "NS", "N", "n_symp","ent","f","n_I")) |>                   
     pivot_wider(names_from = var, values_from = value) |>
     mutate(lhs_num = i) |>                          
     relocate(lhs_num)
    
    # Output balanced partner numbers in 2030
    df_c[[i]] <- result |>
      filter(var %in% c("c","c_ad","c_ad_strat"), year == 2030) |>
      mutate(lhs_num = i) |>                          
      relocate(lhs_num)
    
    # Output per partnership transmission probability in 2030
    df_kappa[[i]] <- result |>
      filter(var %in% c("kappa_p","chi_avg"), year == 2030) |>
      mutate(lhs_num = i) |>                          
      relocate(lhs_num)
  }
  
  print(i)
  
  # return(list(df_result = df_result, df_prev = df_prev))
  return(list(df_prev = df_prev, df_c = df_c, df_kappa = df_kappa))
}
