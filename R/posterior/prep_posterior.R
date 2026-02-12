library(tidyverse)

# Calculate posterior estimate

## DATA ----

lhs_input <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "param")
prior_samples <- data.table::fread("./parameters/prior_samples.csv", check.names = FALSE) |> mutate(lhs_num = row_number())
prior_prev <- data.table::fread("./parameters/prior_prev.csv")

prior_kappa <- data.table::fread("./parameters/prior_kappa.csv")
prior_c <- data.table::fread("./parameters/prior_c.csv")

## FILTERS ----

# Filter FSW prevalence within plausible min and max
fsw_filter <- function(data, min, max) {
  
  data |>
    filter(sex == 2, risk == 3) |>
    group_by(lhs_num, year, sex, risk) |>
    summarise(prev = sum(I)/sum(NS)) |>
    # identify lhs_num that is within prev range
    group_by(lhs_num) |>
    filter(all(prev >= min & prev <= max)) |>
    distinct(lhs_num) |>
    pull(lhs_num)
  
}

# Filter female overall prevalence within plausible min and max
fm_filter <- function(data, min, max) {
  
  data |>
    filter(sex == 2) |>
    group_by(lhs_num, year, sex) |>
    summarise(prev = sum(I)/sum(NS)) |>
    # identify lhs_num that is within prev range
    group_by(lhs_num) |>
    filter(all(prev >= min & prev <= max)) |>
    distinct(lhs_num) |>
    pull(lhs_num)
  
}

# Filter female overall prevalence within plausible min and max
preg_filter <- function(data, min, max) {
  
  data |>
    filter(sex == 2, gest == 2) |>
    group_by(lhs_num, year, sex) |>
    summarise(prev = sum(I)/sum(NS)) |>
    # identify lhs_num that is within prev range
    group_by(lhs_num) |>
    filter(all(prev >= min & prev <= max)) |>
    distinct(lhs_num) |>
    pull(lhs_num)
  
}

# Filter for plausible ratio between male and female prev
# Overall, and not for specific strata
ratio_filter <- function(data, min, max) {
  
  data |>
    group_by(lhs_num, year, sex) |>
    summarise(prev = sum(I)/sum(NS)) |>
    pivot_wider(names_from = sex, values_from = prev) |>
    mutate(ratio = `1`/`2`) |>
    filter(ratio >= min & ratio <= max) |>
    distinct(lhs_num) |>
    pull(lhs_num)
  
}

# Filter for ratio of symptomatic cases
# Overall, and not for specific strata
# Use sexually active population per calculation in DHS
case_filter <- function(data, min, max) {
  
  data |>
    group_by(lhs_num, year, sex) |>
    summarise(prop = sum(n_symp)/sum(NS)) |>
    pivot_wider(names_from = sex, values_from = prop) |>
    mutate(ratio = `1`/`2`) |>
    filter(ratio >= min & ratio <= max) |>
    distinct(lhs_num) |>
    pull(lhs_num)
  
}

# Apply filters 
fsw_nums        <- fsw_filter(prior_prev |> filter(year == 2030), min = 0.02, max = 0.065)
fm_nums         <- fm_filter(prior_prev |> filter(year == 2030), min = 0.01, max = 0.035)
ratio_nums      <- ratio_filter(prior_prev |> filter(year == 2030), min = 0.6, max = 1)
case_nums       <- case_filter(prior_prev |> filter(year == 2030), min = 0.35, max = 0.6)

final_nums <- intersect(fsw_nums, fm_nums) |> intersect(ratio_nums) |> intersect(case_nums) 

write.csv(final_nums, "./posterior/post_lhs_nums.csv", row.names = FALSE)

length(fsw_nums)/nrow(prior_samples)*100
length(fm_nums)/nrow(prior_samples)*100
length(ratio_nums)/nrow(prior_samples)*100
length(case_nums)/nrow(prior_samples)*100
length(final_nums)/nrow(prior_samples)*100

post_prev <- prior_prev|> filter(lhs_num %in% final_nums)
post_samples <- prior_samples |> filter(lhs_num %in% final_nums)
post_c <- prior_c |> filter(lhs_num %in% final_nums)
post_kappa <- prior_kappa |> filter(lhs_num %in% final_nums)

write.csv(post_samples, file = "./posterior/post_samples.csv", row.names = FALSE)
write.csv(post_prev, file = "./posterior/post_prev.csv", row.names = FALSE)
write.csv(post_c, file = "./posterior/post_c.csv", row.names = FALSE)
write.csv(post_kappa, file = "./posterior/post_kappa.csv", row.names = FALSE)

# POSTERIOR OUTCOME QUINTILES ----

# Calculate outcomes for each posterior parameter set 
post_outcomes <- 
  left_join(
    post_prev |>
      filter(year == 2030) |>
      group_by(lhs_num, sex) |>
      summarise(prev = sum(I)/sum(NS)) |>
      mutate(sex = factor(sex, levels = c(1,2), labels = c("prev_m_overall","prev_fm_overall"))) |>
      pivot_wider(names_from = sex, values_from = prev),
    post_prev |>
      filter(year == 2030, sex == 2, risk == 3) |>
      group_by(lhs_num) |>
      summarise(prev_fsw = sum(I)/sum(NS))) |>
  left_join(
    post_prev |>
      filter(year == 2030) |>
      group_by(lhs_num, sex) |>
      summarise(prev = sum(I)/sum(NS)) |>
      pivot_wider(names_from = sex, values_from = prev) |>
      mutate(ratio_mf_prev = `1`/`2`) |>
      select(!c(`1`,`2`))) |> 
  left_join(
    post_prev |>
      filter(year == 2030) |>
      group_by(lhs_num, sex) |>
      summarise(prop = sum(n_symp)/sum(NS)) |>
      pivot_wider(names_from = sex, values_from = prop) |>
      mutate(ratio_mf_cases = `1`/`2`) |>
      select(!c(`1`,`2`))) |>
  ungroup()

# Outcome distribution and correlation
# Strong positive correlation between:
# prev_fm_overall and prev_fsw: when fsw prev increases -> overall fm prevalence increases
# prev_fm_overall and prev_m_overall -> self imposed
# Skewed distributions for most outcomes
post_outcomes |>
  select(!lhs_num) |>
  GGally::ggpairs() 

# DESCRIBE PRIOR AND POSTERIOR ----
order <- c("upsilon-3","upsilon-5","upsilon-4","upsilon-6",
           "q_low-1", "q_med_ratio-1", "q_high_ratio-1",
           
           "c_low-1","c_med_ratio-1","c_med_ratio-2","c_high_ratio-1","c_high_ratio-2","c_age_ratio-1","c_preg_ratio-1",
           
           "chi_low-1","chi_med_ratio-1","chi_high_ratio-1","chi_age_ratio-1",
           "chi_preg_ratio-1","chi_preg_ratio-2","chi_preg_ratio-3",
           
           "e-1",
           "gamma_low-1","gamma_sex_ratio-1","gamma_med_ratio-1","gamma_high_ratio-1","gamma_age_ratio-1","gamma_preg_ratio-1",
           
           "tau_m_low-1", "tau_m_sex_ratio-1",
           "n_l-1",
           "ap-1","ap-2",
           "zeta-1","zeta-2",
           
           "kappa-1","kappa-2",
           "phi-1","phi-2",
           "pi-1","sigma-1"
)

prior_samples |> select(!lhs_num) |>
  summarise(across(everything(), 
                   list(min = min, max = max, mean = mean, sd = sd,
                        median = ~quantile(., probs = 0.5),
                        iqrlwr = ~quantile(., probs = 0.25),
                        iqrupr = ~quantile(., probs = 0.75)))) |>
  pivot_longer(everything(), names_to = c("var", "stat"),
               names_pattern = "(.+)_(.+)") |>
  pivot_wider(names_from = stat, values_from = value) |>
  mutate(var = factor(var, levels = order)) |>
  arrange(var) |>
  print(n=50)

post_samples |> select(!lhs_num) |>
  summarise(across(everything(), 
                   list(min = min, max = max, mean = mean, sd = sd,
                        median = ~quantile(., probs = 0.5),
                        iqrlwr = ~quantile(., probs = 0.25),
                        iqrupr = ~quantile(., probs = 0.75),
                        cilwr = ~quantile(., probs = 0.025),
                        ciupr = ~quantile(., probs = 0.975),
                        n = ~n()))) |>
  pivot_longer(everything(), names_to = c("var", "stat"),
               names_pattern = "(.+)_(.+)") |>
  pivot_wider(names_from = stat, values_from = value) |>
  mutate(var = factor(var, levels = order),
         median_ci = paste0(sprintf("%.2f", median), " (", sprintf("%.2f", cilwr), "-", sprintf("%.2f", ciupr), ")")) |>
  arrange(var) 

# combined table
prior_samples |> select(!lhs_num) |>
  summarise(across(everything(), list(mean = mean, 
                                      cilwr = ~quantile(., probs = 0.025),
                                      ciupr = ~quantile(., probs = 0.975)))) |>
  pivot_longer(everything(), names_to = c("var", "stat"),
               names_pattern = "(.+)_(.+)") |>
  pivot_wider(names_from = stat, values_from = value) |>
  mutate(Prior = case_when(mean < 10 ~ paste0(sprintf("%.2f", mean)," (",sprintf("%.2f", cilwr),"-",sprintf("%.2f", ciupr),")"),
                           TRUE ~ paste0(sprintf("%.1f", mean)," (",sprintf("%.1f", cilwr),"-",sprintf("%.1f", ciupr),")"))) |>
  select(!c(mean, cilwr, ciupr)) |>
  left_join(
    post_samples |> select(!lhs_num) |>
      summarise(across(everything(), 
                       list(mean = mean, 
                            cilwr = ~quantile(., probs = 0.025),
                            ciupr = ~quantile(., probs = 0.975)))) |>
      pivot_longer(everything(), names_to = c("var", "stat"),
                   names_pattern = "(.+)_(.+)") |>
      pivot_wider(names_from = stat, values_from = value) |>
      mutate(Posterior = case_when(mean < 10 ~ paste0(sprintf("%.2f", mean)," (",sprintf("%.2f", cilwr),"-",sprintf("%.2f", ciupr),")"),
                                   TRUE ~ paste0(sprintf("%.1f", mean)," (",sprintf("%.1f", cilwr),"-",sprintf("%.1f", ciupr),")"))) |>
      select(!c(mean, cilwr, ciupr))) |>
  left_join(lhs_input |> select(variable, var_strat, var_latex), 
            by = c("var"="var_strat")) |>
  relocate(var_latex, variable) |>
  mutate(var = factor(var, levels = order)) |>
  arrange(var) |>
  select(!var) |>
  rename(Parameter = var_latex, Definition = variable) |>
  unite("latex", everything(), sep = " & ", remove=FALSE) |>
  mutate(latex = paste0(latex," \\\\")) |>
  write.csv(".././tables/calibration_results.csv", row.names = FALSE)

# DEMOG PARAMS ----

p_levels <- c("m.l.y.n", "m.l.y.p", "m.l.o.n", "m.l.o.p",
              "m.m.y.n", "m.m.y.p", "m.m.o.n", "m.m.o.p",
              "m.h.y.n", "m.h.y.p", "m.h.o.n", "m.h.o.p",
              "f.l.y.n", "f.l.y.p", "f.l.o.n", "f.l.o.p",
              "f.m.y.n", "f.m.y.p", "f.m.o.n", "f.m.o.p",
              "f.h.y.n", "f.h.y.p", "f.h.o.n", "f.h.o.p")

post_prev |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("m","f")),
         risk = factor(risk, levels = c(1,2,3), labels = c("l","m","h")),
         age = factor(age, levels = c(1,2), labels = c("y","o")),
         gest = factor(gest, levels = c(1,2), labels = c("n","p")),
         p = paste0(sex,".",risk,".",age,".",gest),         
         p = factor(p, levels = p_levels)) |> 
  filter(year %in% c(1970,2030)) |>
  group_by(year,sex,p) |>
  summarise(N = mean(N),
            f = mean(f),
            ent = mean(ent)) |>
  pivot_wider(names_from= year, values_from = c(N,f,ent)) |> 
  print(n=100)

# Need to calculate entrants. Actually start at t = 1, but present at t = 0 for simplicity

df_ent <- read.csv("./parameters/entrants.csv")

risk_prop <- post_prev |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("m","f")),
         risk = factor(risk, levels = c(1,2,3), labels = c("l","m","h"))) |>
  filter(year %in% c(2030)) |>
  group_by(year,sex,risk) |>
  summarise(N = mean(N)) |>
  pivot_wider(names_from = risk, values_from = N) |>
  mutate(tot = l+m+h,
         p_l = l/tot, p_m = m/tot, p_h = h/tot)

entries_1970 <- df_ent |>
  filter(ent_t == 1) |>
  pivot_longer(
    cols = c(entm_y, entf_y),
    names_to = "sex",
    values_to = "ent"
  ) |>
  mutate(sex = if_else(sex == "entm_y", "m", "f")) |>
  left_join(risk_prop, by = "sex") |>
  mutate(
    ent_l = ent * p_l,
    ent_m = ent * p_m,
    ent_h = ent * p_h
  ) |>
  # Reshape to match your p format
  pivot_longer(
    cols = c(ent_l, ent_m, ent_h),
    names_to = "risk",
    values_to = "ent_value"
  ) |>
  mutate(
    risk = case_when(
      risk == "ent_l" ~ "l",
      risk == "ent_m" ~ "m",
      risk == "ent_h" ~ "h"
    ),
    p = paste0(sex, ".", risk, ".y.n")  # Only young, non-pregnant get entries
  ) |>
  select(p, ent_value)

# mortality rates
deaths <- read.csv("./parameters/unpopulation_dataportal_20240308173922.csv") |>
  filter(IndicatorName == "Deaths by 5-year age groups and sex") |>
  select(Time, Sex, Age, Value) |>
  rename(year = Time, sex = Sex, age = Age, value = Value)

pop <- read.csv("./parameters/unpopulation_dataportal_20240308173922.csv") |>
  filter(IndicatorName == "Population by 5-year age groups and sex",
         Age %in% c("15-19","20-24","25-29","30-34","35-39","40-44","45-59")) |>
  select(Time, Sex, Age, Value) |>
  rename(year = Time, sex = Sex, age = Age, value = Value)

mortality <- left_join(
  deaths |>
    mutate(age = case_when(age %in% c("15-19","20-24") ~ "y",
                           age %in% c("25-29","30-34","35-39","40-44","45-59") ~ "o")) |>
    filter(!is.na(age), year %in% c(1970,2030)) |>
    group_by(year, sex, age) |>
    summarise(deaths = sum(value)),
  
  pop |>
    mutate(age = case_when(age %in% c("15-19","20-24") ~ "y",
                           age %in% c("25-29","30-34","35-39","40-44","45-59") ~ "o")) |>
    filter(!is.na(age), year %in% c(1970,2030)) |>
    group_by(year, sex, age) |>
    summarise(pop = sum(value))
) |>
  mutate(mortality = deaths/pop) |>
  select(!c(deaths, pop)) |>
  pivot_wider(names_from = year, values_from = mortality) |>
  mutate(sex = case_when(sex == "Female" ~ "f",sex == "Male" ~ "m")) |>
  rename(mu_1970 = `1970`, mu_2030 = `2030`)

# final table to export
post_prev |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("m","f")),
         risk = factor(risk, levels = c(1,2,3), labels = c("l","m","h")),
         age = factor(age, levels = c(1,2), labels = c("y","o")),
         gest = factor(gest, levels = c(1,2), labels = c("n","p")),
         p = paste0(sex,".",risk,".",age,".",gest),         
         p = factor(p, levels = p_levels)) |> 
  filter(year %in% c(1970,2030)) |>
  group_by(year,sex,gest,age,p) |>
  summarise(N = mean(N),
            f = mean(f),
            ent = mean(ent)) |>
  ungroup() |>
  mutate(f = case_when(gest == "p" ~ 0, TRUE~f)) |>
  pivot_wider(names_from= year, values_from = c(N,f,ent)) |>
  left_join(entries_1970, by = "p") |>
  mutate(ent_1970 = coalesce(ent_value, ent_1970)) |>
  left_join(mortality, by = c("sex","age")) |>
  select(!c(ent_value, sex, gest,age)) |>
  mutate(p = factor(p, levels = p_levels)) |>
  arrange(p) |>
  mutate(
    across(starts_with("N_"), ~sprintf("%.0f", .)),
    across(starts_with("f_"), ~sprintf("%.2f", .)),
    across(starts_with("ent_"), ~sprintf("%.0f", .)),
    across(starts_with("mu_"), ~sprintf("%.4f", .))
  ) |>
  mutate(latex = paste0(p, " & \\num{", N_1970, "} & \\num{", N_2030,
                        "} & \\num{", ent_1970, "} & \\num{", ent_2030,
                        "} & \\num{", f_1970, "} & \\num{", f_2030,
                        "} & \\num{", mu_1970, "} & \\num{", mu_2030,
                        "} \\\\")) |>
  write.csv(".././tables/param_demog.csv", row.names = FALSE)

