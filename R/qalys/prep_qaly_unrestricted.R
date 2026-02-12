# Calculate QALY losses for scenario with unrestricted POCT availability

library(tidyverse)
library(patchwork)
library(data.table)

# IMPORT DATA ----

# Model runs to keep
lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

# Number cases and person-years per scenario and model run
df_intervention <- data.table::fread("./posterior/post_intervention_max.csv") |>
  filter(lhs_num %in% lhsnum)

# QALY input parameters
qaly <- readxl::read_excel("./qalys/qalys.xlsx", sheet = "qaly") |>
  select(param_full, mean) |>
  pivot_wider(names_from = param_full, values_from = mean)

# UN WPP 2024 projections for 2030
demog <- read.csv("./parameters/unpopulation_dataportal_20250204122450.csv") |>
  distinct() |>
  rename(var = IndicatorShortName, year = Time, sex = Sex, age_group = Age, 
         age_start = AgeStart, age_end = AgeEnd, value = Value) |>
  select(var, year, sex, age_group, age_start, age_end, value) |>
  mutate(var_short = case_when(
    var == "Age-specific mortality rates by age groups and by sex" ~ "mortality_rate",
    var == "Annual population by 5-year age groups and by sex" ~ "population",
    var == "Deaths by age and sex - abridged" ~ "deaths",
    var == "Life expectancy by sex (at birth)" ~ "life_expectancy",
    # Probability of an individual age x dying before the end of the age interval (x, x+n) where n is the length of the interval
    var == "Probability of dying by age groups and by sex" ~ "prob_death")) |>
  filter(!is.na(var_short))

# PREP INFECTION DATA ----

df_infect <- df_intervention[
  year >= 2024,
  lapply(.SD, sum),
  by = .(lhs_num, year, test_strat, test_num, test_pop, sex, gest, age, risk),
  .SDcols = !c("sex2", "risk2", "age2", "gest2")] |>
  mutate(test_strat = factor(test_strat, levels = c("h","y"),
                             labels = c("Screening", "Diagnostic testing")),
         sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
         gest = factor(gest, levels = c(1,2), labels = c("Non-pregnant", "Pregnant")),
         age = factor(age, levels = c(1,2), labels = c("Young", "Old")),
         risk = factor(risk, levels = c(1,2,3), labels = c("Low", "Medium", "High"))) |>
  mutate(
    # number incident infections in each pathway
    case_symp_treat = n_My + n_Ty + n_Th_symp,
    case_symp_untreat = n_Y + n_Z - case_symp_treat,
    case_asymp_treat = n_Th_asymp,
    case_asymp_untreat = n_X - case_asymp_treat,
    # person-years in each pathway (equivalent to incident cases × duration)
    # note difference between person-years and incidence:
    # when someone moves from Y to other compartments:
    # person-years (py_Y) only accumulate for the time they spent in Y
    # once they move to My/Ty/Ry/Th/Rh, they stop contributing to py_Y and start contributing to py_My/py_Ty/etc.
    py_symp_treat = py_My + py_Ty + py_Th_symp,
    py_symp_untreat = py_Y + py_Z + py_Ry + py_Rh_symp,
    py_asymp_treat = py_Th_asymp,
    py_asymp_untreat = py_X + py_Rh_asymp,
    # tests used
    test_used = case_when(test_strat == "Screening" ~ n_test_h,
                          test_strat == "Diagnostic testing" ~ n_test_y + n_test_o)) |>
  filter(!(sex == "Male" & gest == "Pregnant"))

# Calculate annual accumulation for py
df_infect <- df_infect |>
  group_by(lhs_num, test_pop, test_strat, test_num, sex, gest, age, risk) |>
  mutate(
    annual_py_symp_treat = py_symp_treat - lag(py_symp_treat),
    annual_py_symp_untreat = py_symp_untreat - lag(py_symp_untreat),
    annual_py_asymp_treat = py_asymp_treat - lag(py_asymp_treat),
    annual_py_asymp_untreat = py_asymp_untreat - lag(py_asymp_untreat)
  ) |>
  filter(year != 2024) |>  # Now remove 2024 after calculating differences
  ungroup()

# pivot longer to column for cases and person-years
df <- df_infect |>
  select(lhs_num, test_pop, test_strat, test_num, year, sex, gest, age, risk,
         starts_with("case_symp"), starts_with("case_asymp"),
         starts_with("py_symp"), starts_with("py_asymp"),
         starts_with("annual_py_")) |>
  pivot_longer(cols = starts_with(c("case_", "py_", "annual_py_")),
               names_to = c(".value", "pathway"),
               names_pattern = "(case|py|annual_py)_(.+)",
               values_to = c("cases", "py","annual_py")) 

# Split pathway column using DT since separate() took too long
setDT(df)[, c("inf", "tx") := tstrsplit(pathway, "_", fixed=TRUE)]

glimpse(df)

# test_used is separate because pivot longer above repeats the data across tx and symp status (confusing)
test_used <- df_infect |>
  select(lhs_num, test_pop, test_strat, test_num, test_used, year, sex, gest, age, risk) 

# remove df_intervention as very large
rm(df_intervention)

# CALC YEARS LIVED ----
# For conditions > 1 year, calculated expected years lived with condition
# Accounts for death during this period
# Discounts by 3% to present year

## Chronic conditions ----

## YOUNG GROUP
# Gets infected at mid-point of 15 - 24 year age group (age 19.5 years-> 20)
# Gets condition 5 years later (age 24.5 years -> 25)
# Condition lasts 10 years (age 25 - 34)
yl_young <- demog |>
  select(!var) |>
  filter(var_short %in% c("prob_death"), 
         sex == "Female",
         age_group %in% c("20-24","25-29","30-34")) |>
  pivot_wider(names_from = var_short, values_from = value) |>
  mutate(years_20_tomid = c(2.5, 7.5, 12.5),
         multiplier = 1 - prob_death,
         # Use cumprod() for cumulative multiplication
         # prob_survive = prob_survive prev age group * (1 - prob_death prev age group)
         prob_survive = c(1, cumprod(multiplier[-n()])),
         # years lived = 5 × (probability survive to start of interval) × (1 - probability of dying in interval/2)
         # assumes that deaths occur on average halfway through each interval
         years_lived_undiscount = prob_survive * (prob_death * 2.5 + (1 - prob_death) * 5),
         # discount factor = 1/(1.03^(years from age 20 to midpoint))
         years_lived_discount = years_lived_undiscount * 1 / (qaly$r^years_20_tomid)) 

qaly$yl_fy_chronic_undiscount = yl_young$years_lived_undiscount[2] + yl_young$years_lived_undiscount[3]
qaly$yl_fy_chronic_discount = yl_young$years_lived_discount[2] + yl_young$years_lived_discount[3]

## OLD GROUP
# Gets infected at mid-point of 25 - 49 year age group (age 37 years)
# Gets condition 5 years later (age 42 years)
# Condition lasts 10 years (age 42 - 51)
yl_old <- demog |>
  select(!var) |>
  filter(var_short %in% c("prob_death"), 
         sex == "Female",
         age_group %in% c("35-39","40-44","45-49","50-54")) |>
  pivot_wider(names_from = var_short, values_from = value) |>
  mutate(years_37_tomid = c(0.5, 5.5, 10.5, 15.5),
         multiplier = 1 - prob_death,
         # Use cumprod() for cumulative multiplication
         prob_survive = c(1, cumprod(multiplier[-n()])),
         years_lived_undiscount = prob_survive * (prob_death * 2.5 + (1 - prob_death) * 5),
         years_lived_discount = years_lived_undiscount * 1 / (qaly$r^years_37_tomid)) 

qaly$yl_fo_chronic_undiscount = yl_old$years_lived_undiscount[2]*3/5 + yl_old$years_lived_undiscount[3] + yl_old$years_lived_undiscount[4]*2/5
qaly$yl_fo_chronic_discount = yl_old$years_lived_discount[2]*3/5 + yl_old$years_lived_discount[3] + yl_old$years_lived_discount[4]*2/5

## Pregnant ----

## MOTHER YOUNG
# Adverse pregnancy event happens at mid-point of 15 - 24 year age group (age 19.5 years-> 20)
# "Condition" lasts 10 years (age 20 - 30)
yl_young <- demog |>
  select(!var) |>
  filter(var_short %in% c("prob_death"), 
         sex == "Female",
         age_group %in% c("20-24","25-29","30-34")) |>
  pivot_wider(names_from = var_short, values_from = value) |>
  mutate(years_20_tomid = c(2.5, 7.5, 12.5),
         multiplier = 1 - prob_death,
         # Use cumprod() for cumulative multiplication
         prob_survive = c(1, cumprod(multiplier[-n()])),
         years_lived_undiscount = prob_survive * (prob_death * 2.5 + (1 - prob_death) * 5),
         years_lived_discount = years_lived_undiscount * 1 / (qaly$r^years_20_tomid)) 

qaly$yl_fy_preg_undiscount = yl_young$years_lived_undiscount[1] + yl_young$years_lived_undiscount[2]
qaly$yl_fy_preg_discount = yl_young$years_lived_discount[1] + yl_young$years_lived_discount[2]

## MOTHER OLD
# Adverse pregnancy event happens at mid-point of 25 - 49 year age group (age 37 years)
# "Condition" lasts 10 years (age 37 - 47)
yl_old <- demog |>
  select(!var) |>
  filter(var_short %in% c("prob_death"), 
         sex == "Female",
         age_group %in% c("35-39","40-44","45-49")) |>
  pivot_wider(names_from = var_short, values_from = value) |>
  mutate(years_37_tomid = c(0, 5, 10),
         multiplier = 1 - prob_death,
         # Use cumprod() for cumulative multiplication
         prob_survive = c(1, cumprod(multiplier[-n()])),
         years_lived_undiscount = prob_survive * (prob_death * 2.5 + (1 - prob_death) * 5),
         years_lived_discount = years_lived_undiscount * 1 / (qaly$r^years_37_tomid)) 

qaly$yl_fo_preg_undiscount = yl_old$years_lived_undiscount[1]*3/5 + yl_old$years_lived_undiscount[2] + yl_old$years_lived_undiscount[3]*2/5
qaly$yl_fo_preg_discount = yl_old$years_lived_discount[1]*3/5 + yl_old$years_lived_discount[2] +  yl_old$years_lived_discount[3]*2/5

## Newborn ----
# Years lived until life expectancy
# 64.8766 -> 65 years on average
# demog |>
#   select(!var) |>
#   filter(var_short %in% c("life_expectancy")) 

yl_n <- demog |>
  select(!var) |>
  filter(var_short %in% "prob_death",
         sex == "Both sexes",
         age_start < 65) |>
  pivot_wider(names_from = var_short, values_from = value) |>
  mutate(years_0_tomid = c(0.5, 3.5, seq(7.5, 62.5, by = 5)),
         multiplier = 1 - prob_death,
         # Use cumprod() for cumulative multiplication
         prob_survive = c(1, cumprod(multiplier[-n()])),
         years_lived_undiscount = case_when(age_group == "0" ~ 1 * prob_survive * (1 - prob_death/2),
                                            age_group == "1-4" ~ 4 * prob_survive * (1 - prob_death/2),
                                            TRUE ~ prob_survive * (prob_death * 2.5 + (1 - prob_death) * 5)),
         years_lived_discount = years_lived_undiscount * 1 / (qaly$r^years_0_tomid)) 

qaly$yl_n_life_undiscount <- yl_n  |> group_by() |> summarise(n = sum(years_lived_undiscount)) |> pull(n)
qaly$yl_n_life_discount <- yl_n  |> group_by() |> summarise(n = sum(years_lived_discount)) |> pull(n)

qaly$yl_n_10_undiscount <- yl_n |> filter(age_start < 10) |> group_by() |> summarise(n = sum(years_lived_undiscount)) |> pull(n)
qaly$yl_n_10_discount <- yl_n |> filter(age_start < 10) |> group_by() |> summarise(n = sum(years_lived_discount)) |> pull(n)

glimpse(qaly)

# CALC QALYS ----
# QALYs lost = number infections × prob sequela × duration of sequela × (1 - utility of sequela)

# For urethritis, ng, pid and pid sequelae:
# originally used incident cases * duration of infection to determine person-years
# but can rather get person-years directly from the model without needing complex duration calcs

## Men ----
df_qaly_male <- df |>
  filter(sex == "Male") |>
  mutate(
    qaly_m_urethritis = case_when(
      inf == "symp" ~ annual_py * qaly$prob_urth_symp * (1 - qaly$util_ng_m_symp),
      inf == "asymp" ~ annual_py * qaly$prob_urth_asymp * (1 - qaly$util_ng_asymp),
      TRUE ~ 0),
    qaly_m_eds = case_when(
      # note that EDS calc is independent of duration of infection (unlike other sequalae - see paper)
      tx == "treat" ~ case * qaly$prob_eds_treat * qaly$dur_eds * (1 - qaly$util_eds),
      tx == "untreat" ~ case * qaly$prob_eds_untreat * qaly$dur_eds * (1 - qaly$util_eds),
      TRUE ~ 0)
  )

## Non-pregnant women ----

df_qaly_female_nonpreg <- df |>
  filter(sex == "Female" & gest == "Non-pregnant") |>
  mutate(
    qaly_f_ng = case_when(
      inf == "symp" ~ annual_py * (1 - qaly$util_ng_f_symp),
      inf == "asymp" ~ annual_py * (1 - qaly$util_ng_asymp),
      TRUE ~ 0),
    
    # PID
    # note that there are two duration terms here
    # prob PID  = annual prob PID * duration of infection
    # qaly = case * prob PID * duration PID 
    # not differentiated for symp/asymp or treat/untreat
    qaly_f_pid = annual_py * qaly$pid_annual * qaly$dur_pid * (1 - qaly$util_pid),
    
    qaly_f_ep = annual_py * qaly$pid_annual * qaly$prob_ep * qaly$dur_ep * (1 - qaly$util_ep),
    
    qaly_f_cpp = case_when(
      age == "Young" ~ 
        # number of CPP cases
        annual_py * qaly$pid_annual * qaly$prob_cpp * 
        # discounted years with CPP 
        qaly$yl_fy_chronic_discount *
        # disutility
        (1 - qaly$util_cpp),
      age == "Old" ~ 
        # number of CPP cases
        annual_py * qaly$pid_annual * qaly$prob_cpp * 
        # discounted years with CPP 
        qaly$yl_fo_chronic_discount *
        # disutility
        (1 - qaly$util_cpp),
      TRUE ~ 0),
    
    qaly_f_tfi = case_when(
      # same method as CPP above
      age == "Young" ~ 
        annual_py * qaly$pid_annual * qaly$prob_tfi * qaly$yl_fy_chronic_discount * (1 - qaly$util_tfi),
      age == "Old" ~ 
        annual_py * qaly$pid_annual * qaly$prob_tfi * qaly$yl_fo_chronic_discount * (1 - qaly$util_tfi),
      TRUE ~ 0),
    
    # joint sequelae
    # joint probability is multiplication of sequelae specific probabilities
    # joint utility applies to shortest duration of the two sequelae (joint util is multiplication of sequelae specific utils)
    # single utility then applies to remaining duration for longer sequela
    # e.g. QALY lost for CPP and EP = (1 - utility_c * utility_e) * dur_e  +  (1 - util_c) * (dur_c - dur_e)
    # see 2015 paper for these details
    
    qaly_f_cpp_tfi = case_when(
      # both chronic conditions with same duration, so apply joint utility for entire duration
      age == "Young"  ~ 
        annual_py * qaly$pid_annual * qaly$prob_cpp * qaly$prob_tfi * qaly$yl_fy_chronic_discount * (1 - qaly$util_cpp * qaly$util_tfi),
      age == "Old"  ~ 
        annual_py * qaly$pid_annual * qaly$prob_cpp * qaly$prob_tfi * qaly$yl_fo_chronic_discount * (1 - qaly$util_cpp * qaly$util_tfi),
      TRUE ~ 0),
    
    qaly_f_cpp_ep = case_when(
      age == "Young" ~ 
        # Part 1: joint utility loss during EP duration (shorter) + 
        # Part 2: CPP-only utility loss for the remaining chronic duration
        annual_py * qaly$pid_annual * qaly$prob_cpp * qaly$prob_ep * 
        (qaly$dur_ep * (1 - qaly$util_cpp * qaly$util_ep) + (qaly$yl_fy_chronic_discount - qaly$dur_ep) * (1 - qaly$util_cpp)),
      age == "Old" ~ 
        annual_py * qaly$pid_annual * qaly$prob_cpp * qaly$prob_ep * 
        (qaly$dur_ep * (1 - qaly$util_cpp * qaly$util_ep) + (qaly$yl_fo_chronic_discount - qaly$dur_ep) * (1 - qaly$util_cpp)),
      TRUE ~ 0),
    
    qaly_f_ep_tfi = case_when(
      age == "Young" ~  
        # Part 1: joint utility loss during EP duration (shorter)
        # Part 2: TFI-only utility loss for the remaining chronic duration
        annual_py * qaly$pid_annual * qaly$prob_ep * qaly$prob_tfi *
        (qaly$dur_ep * (1 - qaly$util_ep * qaly$util_tfi) + (qaly$yl_fy_chronic_discount - qaly$dur_ep) * (1 - qaly$util_tfi)),
      age == "Old" ~  
        annual_py * qaly$pid_annual * qaly$prob_ep * qaly$prob_tfi *
        (qaly$dur_ep * (1 - qaly$util_ep * qaly$util_tfi) + (qaly$yl_fo_chronic_discount - qaly$dur_ep) * (1 - qaly$util_tfi)),
      TRUE ~ 0),
    
    qaly_f_cpp_tfi_ep = case_when(
      age == "Young" ~ 
        # Part 1: three-way joint utility loss during EP duration (shortest)
        # Part 2: CPP+TFI joint utility loss for the remaining chronic duration
        annual_py * qaly$pid_annual * qaly$prob_cpp * qaly$prob_tfi * qaly$prob_ep * 
        (qaly$dur_ep * (1 - qaly$util_cpp * qaly$util_tfi * qaly$util_ep) + (qaly$yl_fy_chronic_discount - qaly$dur_ep) * (1 - qaly$util_cpp * qaly$util_tfi)),
      age == "Old" ~ 
        annual_py * qaly$pid_annual * qaly$prob_cpp * qaly$prob_tfi * qaly$prob_ep * 
        (qaly$dur_ep * (1 - qaly$util_cpp * qaly$util_tfi * qaly$util_ep) + (qaly$yl_fo_chronic_discount - qaly$dur_ep) * (1 - qaly$util_cpp * qaly$util_tfi)),
      TRUE ~ 0))

## Pregnant women ----
# Note that utilities use incident cases and not person-years
# Probability of event occuring is once per infection, not for the duration of infection
# i.e if infection lasted 2 years wouldn't have twice the risk of stillbirth (according to how probability was calculated)

# df_qaly_female_preg <- df |>
#   filter(sex == "Female" & gest == "Pregnant") |>
#   mutate(
#     qaly_f_ng = case_when(
#       inf == "symp" ~ annual_py * (1 - qaly$util_ng_f_symp),
#       inf == "asymp" ~ annual_py * (1 - qaly$util_ng_asymp),
#       TRUE ~ 0),
#     qaly_f_still = case_when(
#       tx == "untreat" & age == "Young" ~ annual_py * qaly$prob_still * qaly$yl_fy_preg_discount * (1 - qaly$util_f_still),
#       tx == "untreat" & age == "Old" ~ annual_py * qaly$prob_still * qaly$yl_fo_preg_discount * (1 - qaly$util_f_still)),
#     qaly_f_lbw = case_when(
#       tx == "untreat" & age == "Young" ~ annual_py * (1 - qaly$prob_still) * qaly$prob_lbw * qaly$yl_fy_preg_discount * (1 - qaly$util_f_lbw),
#       tx == "untreat" & age == "Old" ~ annual_py * (1 - qaly$prob_still) * qaly$prob_lbw * qaly$yl_fo_preg_discount * (1 - qaly$util_f_lbw)),
#     qaly_f_conj = case_when(
#       tx == "untreat" ~ annual_py * (1 - qaly$prob_still) * qaly$prob_conj * qaly$dur_f_conj * (1 - qaly$util_f_conj)),
#     qaly_f_pneu = case_when(
#       tx == "untreat" ~ annual_py * (1 - qaly$prob_still) * qaly$prob_pneu * qaly$dur_f_pneu * (1 - qaly$util_f_pneu))
#     )

df_qaly_female_preg <- df |>
  filter(sex == "Female" & gest == "Pregnant") |>
  mutate(
    qaly_f_ng = case_when(
      inf == "symp" ~ annual_py * (1 - qaly$util_ng_f_symp),
      inf == "asymp" ~ annual_py * (1 - qaly$util_ng_asymp),
      TRUE ~ 0),
    qaly_f_still = case_when(
      tx == "untreat" & age == "Young" ~ case * qaly$prob_still * qaly$yl_fy_preg_discount * (1 - qaly$util_f_still),
      tx == "untreat" & age == "Old" ~ case * qaly$prob_still * qaly$yl_fo_preg_discount * (1 - qaly$util_f_still)),
    qaly_f_lbw = case_when(
      tx == "untreat" & age == "Young" ~ case * (1 - qaly$prob_still) * qaly$prob_lbw * qaly$yl_fy_preg_discount * (1 - qaly$util_f_lbw),
      tx == "untreat" & age == "Old" ~ case * (1 - qaly$prob_still) * qaly$prob_lbw * qaly$yl_fo_preg_discount * (1 - qaly$util_f_lbw)),
    qaly_f_conj = case_when(
      tx == "untreat" ~ case * (1 - qaly$prob_still) * qaly$prob_conj * qaly$dur_f_conj * (1 - qaly$util_f_conj)),
    qaly_f_pneu = case_when(
      tx == "untreat" ~ case * (1 - qaly$prob_still) * qaly$prob_pneu * qaly$dur_f_pneu * (1 - qaly$util_f_pneu))
  )


## Newborn ----
# Note that utilities use incident cases and not person-years

# df_qaly_newborn <- df |>
#   filter(sex == "Female" & gest == "Pregnant", tx == "untreat") |>
#   mutate(
#     qaly_n_still = annual_py * qaly$prob_still * qaly$yl_n_life_discount * (1 - qaly$util_n_still),
#     # remaining sequelae conditional on live birth
#     qaly_n_lbw = annual_py * (1 - qaly$prob_still) * qaly$prob_lbw * qaly$yl_n_10_discount * (1 - qaly$util_n_lbw),
#     qaly_n_conj = annual_py * (1 - qaly$prob_still) * qaly$prob_conj * qaly$dur_n_conj * (1 - qaly$util_n_conj),
#     qaly_n_pneu = annual_py * (1 - qaly$prob_still) * qaly$prob_pneu * qaly$dur_n_pneu * (1 - qaly$util_n_pneu)
#     )
df_qaly_newborn <- df |>
  filter(sex == "Female" & gest == "Pregnant", tx == "untreat") |>
  mutate(
    qaly_n_still = case * qaly$prob_still * qaly$yl_n_life_discount * (1 - qaly$util_n_still),
    # remaining sequelae conditional on live birth
    qaly_n_lbw = case * (1 - qaly$prob_still) * qaly$prob_lbw * qaly$yl_n_life_discount * (1 - qaly$util_n_lbw),
    qaly_n_conj = case * (1 - qaly$prob_still) * qaly$prob_conj * qaly$dur_n_conj * (1 - qaly$util_n_conj),
    qaly_n_pneu = case * (1 - qaly$prob_still) * qaly$prob_pneu * qaly$dur_n_pneu * (1 - qaly$util_n_pneu)
  )

# FORMAT AND SAVE ----
df_qaly <- 
  full_join(df_qaly_female_nonpreg,
            df_qaly_female_preg) |>
  full_join(df_qaly_male) |>
  full_join(df_qaly_newborn)

# Summarise overall 
df_qaly_tot <- df_qaly |>
  pivot_longer(cols = starts_with("qaly"), names_to = "var", values_to = "qaly") |>
  filter(!qaly == 0, !is.na(qaly)) |>
  group_by(lhs_num, test_pop, test_strat, test_num) |>
  summarise(qaly = sum(qaly)) |>
  ungroup() |>
  group_by(lhs_num, test_pop, test_strat) |>
  mutate(qaly_baseline = qaly[test_num == 0],
         qaly_gained = qaly_baseline - qaly) |>
  ungroup() |>
  filter(!test_num == 0) |>
  select(!qaly) |>
  left_join(df_qaly |>
              group_by(lhs_num, test_pop, test_strat, test_num) |>
              summarise(case = sum(case),
                        annual_py = sum(annual_py))) |>
  left_join(test_used |>
              group_by(lhs_num, test_pop, test_strat, test_num) |>
              summarise(test_used = sum(test_used))) |>
  select(lhs_num, test_pop, test_strat, test_num, test_used, case, annual_py, qaly_baseline, qaly_gained)


data.table::fwrite(df_qaly_tot, "./qalys/impact_qaly_unrestricted_tot.csv") 
