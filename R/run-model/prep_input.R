## PREPARE MODEL INPUTS

library(tidyverse)

# 1. Input parameters ----
# If type == "baseline", process df_excel to prepare parameters and calculate arrays
# If type == "update", process df_input with new param values to calculate arrays 
prep_param <- function(df, type) {
  
  if (type == "baseline") {
    
    df <- df |>
      mutate(sex = factor(sex, levels = c("male","female"), labels = c(1,2)),
             risk = factor(risk, levels = c("low","med","high"), labels = c(1,2,3)),
             age = factor(age, levels = c("young","old"), labels = c(1,2)),
             gest = factor(gest, levels = c("non-pregnant","pregnant"), labels = c(1,2))) |>
      select(!c(category, variable, unit, notes, min, max, location_array, var_strat, 
                var_ggplot, var_latex, var_ggplot_old, var_latex_old))
    
    # Split dataframe by variable
    variable_list <- split(df, df$var)
    
    # Create a list to store the arrays
    df_list <- list()
    
    # Iterate over each variable
    for (var in names(variable_list)) {
      
      # Extract the dataframe for the current variable
      dat <- variable_list[[var]]
      
      # Remove empty columns
      dat <- dat[colSums(is.na(dat)) == 0]
      
      # Get dimensions
      dim <- case_when(ncol(dat) == 6 & "risk" %in% names(dat) ~ "c(2,3,2,2)",
                       ncol(dat) == 5 & "risk" %in% names(dat) ~ "c(2,3,2)",
                       ncol(dat) == 5 & !"risk" %in% names(dat) ~ "c(2,2,2)",
                       ncol(dat) == 4 & "risk" %in% names(dat) ~ "c(2,3)",
                       ncol(dat) == 4 & !"risk" %in% names(dat) ~ "c(2,2)",
                       ncol(dat) == 3 ~ NA,     # vector
                       ncol(dat) == 1 ~ NA)     # value
      
      if(!is.na(dim)){
        
        # Create array
        array <- array(
          xtabs(value ~ ., data = dat),
          dim = eval(parse(text = dim)))
        
        # Assign array to list
        df_list[[var]] <- array
        
      } else if (is.na(dim)) {
        
        # Assign vector
        df_list[[var]] <- dat$value
        
      } 
      
    }
    
  } else if (type == "update") {
    
    # If updating, df is already structured as df_list
    df_list <- df
    
  } 
  
  ## Now process or update variables
  
  # q: probability of being sexually active
  # calculated for young non pregnant, and set as 1 for old and pregnant
  q <- df_list$q_low
  q[2,1,1,1] <- q[1,1,1,1] * df_list$q_sex_ratio
  q[,2,1,1] <- q[,1,1,1] * df_list$q_med_ratio
  q[,3,1,1] <- q[,2,1,1] * df_list$q_high_ratio
  q[,,2,1] <- df_list$q_old
  q[2,,,2] <- df_list$q_preg
  df_list$q <- q
  
  # Upsilon: calculate without ratio
  upsilon <- df_list$upsilon
  upsilon[,1] <- 1 - upsilon[,2] - upsilon[,3]
  df_list$upsilon <- upsilon
  
  # c: Number sex partners
  c <- df_list$c_low
  c[2,1,1,1] <- c[1,1,1,1] * df_list$c_sex_ratio
  c[,1,2,1] <- c[,1,1,1] * df_list$c_age_ratio[] # older low risk 
  c[2,1,,2] <- c[2,1,,1] * df_list$c_preg_ratio[] # pregnant low risk
  c[,2,,] <- c[,1,,] * df_list$c_med_ratio[] # medium risk for all strata
  c[,3,,] <- c[,2,,] * df_list$c_high_ratio[] # high risk for all strata
  df_list$c <- c
  
  # chi: Condom use
  chi <- df_list$chi_low
  chi[2,1,1,1] <- chi[1,1,1,1] 
  chi[,1,2,1] <- chi[,1,1,1] * df_list$chi_age_ratio[] # older low risk 
  chi[,2,,] <- chi[,1,,] * df_list$chi_med_ratio[] # medium risk, all age, non-preg
  chi[,3,,] <- chi[,2,,] * df_list$chi_high_ratio[] # high risk, all age, non-preg
  chi[2,,,2] <- chi[2,,,1] * df_list$chi_preg_ratio[] # pregnant, by risk group
  df_list$chi <- chi
  
  # hc: Proportion accessing PHC or ANC
  hc <- df_list$hc_low
  hc[,1,2,1] <- hc[,1,1,1] * df_list$hc_age_ratio[] # older low risk 
  hc[2,1,,2] <- hc[2,1,,1] * df_list$hc_preg_ratio[] # pregnant low risk
  hc[,2,,] <- hc[,1,,] * df_list$hc_med_ratio[] # medium risk for all strata
  hc[,3,,] <- hc[,1,,] * df_list$hc_high_ratio[] # high risk for all strata
  df_list$hc <- hc
  
  # gamma: Proportion accessing STI services
  gamma <- df_list$gamma_low
  gamma[2,1,1,1] <- gamma[1,1,1,1] * df_list$gamma_sex_ratio
  gamma[,1,2,1] <- gamma[,1,1,1] * df_list$gamma_age_ratio # older low risk 
  gamma[,2,,] <- gamma[,1,,] * df_list$gamma_med_ratio # medium risk, all age, non-preg
  gamma[,3,,] <- gamma[,2,,] * df_list$gamma_high_ratio # high risk, all age, non-preg
  gamma[2,,,2] <- gamma[2,,,1] * df_list$gamma_preg_ratio[] # pregnant, by risk group
  df_list$gamma <- gamma
  
  # tau_m: Rate of seeking treatment for symptomatic infection
  tau_m <- df_list$tau_m_low
  tau_m[2,1,1,1] <- tau_m[1,1,1,1] * df_list$tau_m_sex_ratio
  tau_m[,2,,] <- tau_m[,1,,] * df_list$tau_m_med_ratio[] # medium risk for all strata
  tau_m[1,3,,] <- tau_m[1,1,,] * df_list$tau_m_high_ratio[] # high risk for all strata
  tau_m[2,3,,] <- tau_m[1,1,,] 
  tau_m[,,2,] <- tau_m[,,1,] * df_list$tau_m_age_ratio[] # older low risk 
  tau_m[2,,,2] <- tau_m[2,,,1] * df_list$tau_m_preg_ratio[] # pregnant low risk
  df_list$tau_m <- tau_m
  
  # pi and sigma duplicated per sex
  df_list$pi <- rep(df_list$pi,2)
  df_list$sigma <- rep(df_list$sigma,2)
  
  # number sex acts
  df_list$n_m <- df_list$n_l * df_list$n_m_ratio
  df_list$n_h <- df_list$n_m * df_list$n_h_ratio
  
  # VARIABLE TO CHECK FOR ANY PROBABILITIES > 1
  # df_list$error_check <- if ((any(q > 1) | any(upsilon > 1 | upsilon < 0) | any(chi > 1) | 
  #                     any(hc > 1) | any(gamma > 1))) TRUE else FALSE
  
  df_list$error_check <- c()
  
  if (any(q > 1)) {
    df_list$error_check <- c(df_list$error_check, "q")
  }
  
  if (any(upsilon > 1 | upsilon < 0)) {
    df_list$error_check <- c(df_list$error_check, "upsilon")
  }
  
  if (any(chi > 1)) {
    df_list$error_check <- c(df_list$error_check, "chi")
  }
  
  if (any(hc > 1)) {
    df_list$error_check <- c(df_list$error_check, "hc")
  }
  
  if (any(gamma > 1)) {
    df_list$error_check <- c(df_list$error_check, "gamma")
  }
  
  # If no errors were found, set error_check to NA
  if (length(df_list$error_check) == 0) {
    df_list$error_check <- NA_character_
  }
  
  return(df_list)
  
}

# 2. Initial conditions ----
# Use input parameters (upsilon and q) to divide N by risk group and calculate S
# Assign initial values to other state variables 
prep_init <- function(df_excel_init, df_param) {
  
  init <- df_excel_init |>
    filter(var == "N") |>
    mutate(sex = factor(sex, levels = c("male","female"), labels = c(1,2)),
           risk = factor(risk, levels = c("low","med","high"), labels = c(1,2,3)),
           age = factor(age, levels = c("young","old"), labels = c(1,2)),
           gest = factor(gest, levels = c("non-pregnant","pregnant"), labels = c(1,2))) |>
    select(!c(category, variable, unit, notes))
  
  
  # N0
  # N stratified by sex, age, gest
  # Multiply by upsilon: risk group probability by sex
  N0 <- array(xtabs(value ~ ., data = rbind(init |> mutate(risk = 1),
                                            init |> mutate(risk = 2),
                                            init |> mutate(risk = 3))),
       dim = c(2,3,2,2)) 
  
  N0[,1,,] <- N0[,1,,] * df_param$upsilon[,1]
  N0[,2,,] <- N0[,2,,] * df_param$upsilon[,2]
  N0[,3,,] <- N0[,3,,] * df_param$upsilon[,3]
  
  # S0
  # Among young: Multiply N0 by q: proportion sexually active on pop entry
  # Among old: Multiply N0 by 1
  S0 <- N0 * df_param$q

  # Remaining initial parameters
  U0  <- N0 - S0  
  X0  <- N0*0
  Y0  <- N0*0
  Z0  <- N0*0
  
  # initial infections
  X0[1,1,,1] <- N0[1,1,,1]*0.004/2
  Z0[1,1,,1] <- N0[1,1,,1]*0.004/2
  
  X0[1,2,,1] <- N0[1,2,,1]*0.012/2
  Z0[1,2,,1] <- N0[1,2,,1]*0.012/2
  
  X0[1,3,,1] <- N0[1,3,,1]*0.02/2
  Z0[1,3,,1] <- N0[1,3,,1]*0.02/2
  
  X0[2,1,,] <- N0[2,1,,]*0.004/2
  Z0[2,1,,] <- N0[2,1,,]*0.004/2
  
  X0[2,2,,] <- N0[2,2,,]*0.025/2
  Z0[2,2,,] <- N0[2,2,,]*0.025/2
  
  X0[2,3,,] <- N0[2,3,,]*0.04/2
  Z0[2,3,,] <- N0[2,3,,]*0.04/2
  
  S0 <- S0 - X0 - Z0
  
  # check
  # round(S0 + U0 + X0 + Y0 + Z0,0) == round(N0,0)
  
  # add to list
  init_list <- list()
  init_list[["U0"]] <- U0
  init_list[["S0"]] <- S0
  init_list[["X0"]] <- X0
  init_list[["Y0"]] <- Y0
  init_list[["Z0"]] <- Z0
  
  init_list
  
}

