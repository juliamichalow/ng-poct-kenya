prep_df <- function(df) {
  
  df_plot <- df |>
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
    mutate(sex = factor(sex, levels = c(1,2), labels = c("male","female")),
           sex2 = factor(sex2, levels = c(1,2), labels = c("male","female")),
           risk = factor(risk, levels = c(1,2,3), labels = c("low","med","high")),
           risk2 = factor(risk2, levels = c(1,2,3), labels = c("low","med","high")),
           age = factor(age, levels = c(1,2), labels = c("young","old")),
           age2 = factor(age2, levels = c(1,2), labels = c("young","old")),
           gest = factor(gest, levels = c(1,2), labels = c("nonp","preg")),
           gest2 = factor(gest2, levels = c(1,2), labels = c("nonp","preg"))) |>
    select(t, year, var, sex:gest2, value) 
  
  return(df_plot)
  
}


