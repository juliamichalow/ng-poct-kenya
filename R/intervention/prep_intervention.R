# Calculate intervention impact for each scenario

library(tidyverse)
library(data.table)

# DATA PREP AND SAVE ----

post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)

lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

df_intervention_strata <- data.table::fread("./posterior/post_intervention_constrained_strata.csv") |>
  filter(lhs_num %in% lhsnum)

df_intervention_strata |> group_by() |> summarise(n=n_distinct(lhs_num))
df_intervention_strata |> count(test_strat, test_pop, test_num)

# by sex
df_intervention_sex <- df_intervention_strata[year %in% 2024:2030,
                                              lapply(.SD, sum),
                                              by = .(lhs_num, year, test_pop, test_strat, test_num, sex),
                                              .SDcols = !c("risk", "age", "gest", "sex2", "risk2", "age2", "gest2")]

df_intervention_sex <- df_intervention_sex |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
         test_strat = factor(test_strat, levels = c("h","y"),
                             labels = c("Screening", "Diagnostic testing"))) |>
  mutate(eligible = case_when(test_strat == "Screening" ~ n_elig_h,
                              test_strat == "Diagnostic testing" ~ n_elig_y + n_elig_o),
         test_used = case_when(test_strat == "Screening" ~ n_test_h,
                               test_strat == "Diagnostic testing" ~ n_test_y + n_test_o))

fwrite(df_intervention_sex, "./posterior/post_intervention_constrained_sex.csv")

# overall
df_intervention <- df_intervention_strata[year %in% 2024:2030,
                                          lapply(.SD, sum),
                                          by = .(lhs_num, year, test_pop, test_strat, test_num),
                                          .SDcols = !c("sex", "risk", "age", "gest", "sex2", "risk2", "age2", "gest2")]

df_intervention_total <- df_intervention |>
  mutate(test_strat = factor(test_strat, levels = c("h","y"),
                             labels = c("Screening", "Diagnostic testing"))) |>
  mutate(eligible = case_when(test_strat == "Screening" ~ n_elig_h,
                              test_strat == "Diagnostic testing" ~ n_elig_y + n_elig_o),
         test_used = case_when(test_strat == "Screening" ~ n_test_h,
                               test_strat == "Diagnostic testing" ~ n_test_y + n_test_o))

fwrite(df_intervention_total, "./posterior/post_intervention_constrained_total.csv")

# FOR LIMITED NUMBER POCT ----
# IMPACT ----

## Overall----

# Import data if has already been defined and saved above
df_intervention_total <- read.csv("./posterior/post_intervention_constrained_total.csv") 

df_intervention_total |> summarise(n = n_distinct(lhs_num))

# Calculate outcomes:
  # cases averted and prop cases averted relative to baseline, 2025-2030
  # incidence reduction and prop incidence reduction relative to baseline, 2030
df <- left_join(
  df_intervention_total |> 
    filter(year %in% c(2025:2030)) |>
    group_by(lhs_num, test_pop, test_strat, test_num) |>
    summarise(cases = sum(n_I),
              eligible = sum(eligible),
              test_used = sum(test_used)) |>
    ungroup() |>
    group_by(lhs_num, test_pop, test_strat) |>
    # Calculate change relative to baseline (test_num = 0)
    mutate(cases_baseline = cases[test_num == 0],
           cases_averted = cases_baseline - cases,
           prop_cases_averted = cases_averted / cases_baseline,
           test_avail = test_num*6) |>
    ungroup() |>
    filter(test_num != 0) |>
    rename(cases_intervention = cases) |>
    relocate(lhs_num, test_strat, test_pop, test_num, eligible, test_avail, test_used, cases_baseline, cases_intervention, cases_averted, prop_cases_averted),
  
  df_intervention_total |> 
    filter(year == 2030) |>
    mutate(incidence_intervention = n_I/(NS - I)) |>
    group_by(lhs_num, test_pop, test_strat) |>
    # Calculate percentage change relative to baseline (test_num = 0)
    mutate(incidence_baseline = incidence_intervention[test_num == 0],
           prop_incidence_reduction = (incidence_baseline - incidence_intervention) / incidence_baseline) |>
    select(lhs_num, test_strat, test_pop, test_num, incidence_baseline, incidence_intervention, prop_incidence_reduction) |>
    ungroup() |>
    filter(test_num != 0)
)

write.csv(df, "./intervention/impact_intervention_total.csv", row.names = FALSE)

## Stratified ----
p_levels <- c("m.l.y.n", "m.l.o.n", 
              "m.m.y.n", "m.m.o.n",
              "m.h.y.n", "m.h.o.n", 
              "f.l.y.n", "f.l.o.n", "f.l.y.p", "f.l.o.p",
              "f.m.y.n", "f.m.o.n", "f.m.y.p", "f.m.o.p",
              "f.h.y.n", "f.h.o.n", "f.h.y.p", "f.h.o.p")

df_intervention_strata <- data.table::fread("../data/post_intervention_constrained_strata.csv") |>
  filter(year %in% c(2025:2030)) |>
  filter(!(sex ==1 & gest==2)) |>
  mutate(p = interaction(c("m", "f")[sex],
                         c("l", "m", "h")[risk],
                         c("y", "o")[age],
                         c("n", "p")[gest],
                         sep = "."),
         p = factor(p, levels = p_levels),
         test_strat = factor(test_strat, levels = c("h","y"),
                             labels = c("Screening", "Diagnostic testing")),
         eligible = case_when(test_strat == "Screening" ~ n_elig_h,
                            test_strat == "Diagnostic testing" ~ n_elig_y + n_elig_o),
         test_used = case_when(test_strat == "Screening" ~ n_test_h,
                             test_strat == "Diagnostic testing" ~ n_test_y + n_test_o))

df_strat <- left_join(
  df_intervention_strata |> 
    filter(year %in% c(2025:2030)) |>
    group_by(lhs_num, test_pop, test_strat, test_num, p, sex, risk, age, gest) |>
    summarise(cases = sum(n_I),
              eligible = sum(eligible),
              test_used = sum(test_used)) |>
    ungroup() |>
    group_by(lhs_num, test_pop, test_strat, p, sex, risk, age, gest) |>
    # Calculate change relative to baseline (test_num = 0)
    mutate(cases_baseline = cases[test_num == 0],
           cases_averted = cases_baseline - cases,
           prop_cases_averted = cases_averted / cases_baseline,
           test_avail = test_num*6) |>
    ungroup() |>
    filter(test_num != 0) |>
    rename(cases_intervention = cases) |>
    relocate(lhs_num, test_strat, test_pop, test_num, p, sex, risk, age, gest, eligible, test_avail, test_used, cases_baseline, cases_intervention, cases_averted, prop_cases_averted),
  
  df_intervention_strata |> 
    filter(year == 2030) |>
    mutate(incidence_intervention = n_I/(NS - I)) |>
    group_by(lhs_num, test_pop, test_strat, p, sex, risk, age, gest) |>
    # Calculate percentage change relative to baseline (test_num = 0)
    mutate(incidence_baseline = incidence_intervention[test_num == 0],
           prop_incidence_reduction = (incidence_baseline - incidence_intervention) / incidence_baseline) |>
    select(lhs_num, test_strat, test_pop, test_num, p, sex, risk, age, gest, incidence_baseline, incidence_intervention, prop_incidence_reduction) |>
    ungroup() |>
    filter(test_num != 0)
)

fwrite(df_strat, "./intervention/impact_intervention_constrained_strata.csv")

## Annual ----
df_year <- left_join(
  df_intervention_strata |> 
    mutate(test_strat = factor(test_strat, levels = c("h","y"),
                               labels = c("Screening", "Diagnostic testing"))) |>
    mutate(eligible = case_when(test_strat == "Screening" ~ n_elig_h,
                                test_strat == "Diagnostic testing" ~ n_elig_y + n_elig_o),
           test_used = case_when(test_strat == "Screening" ~ n_test_h,
                                 test_strat == "Diagnostic testing" ~ n_test_y + n_test_o)) |>
    filter(year %in% c(2025:2030)) |>
    group_by(lhs_num, test_pop, test_strat, test_num, year) |>
    summarise(cases = sum(n_I),
              eligible = sum(eligible),
              test_used = sum(test_used)) |>
    ungroup() |>
    group_by(lhs_num, test_pop, test_strat, year) |>
    # Calculate change relative to baseline (test_num = 0)
    mutate(cases_baseline = cases[test_num == 0],
           cases_averted = cases_baseline - cases,
           cases_averted = ifelse(cases_averted <0, 0, cases_averted),
           prop_cases_averted = cases_averted / cases_baseline,
           test_avail = test_num) |>
    ungroup() |>
    rename(cases_intervention = cases) |>
    relocate(lhs_num, test_strat, test_pop, test_num, eligible, test_avail, test_used, cases_baseline, cases_intervention, cases_averted, prop_cases_averted),
  
  df_intervention_strata |> 
    mutate(test_strat = factor(test_strat, levels = c("h","y"),
                               labels = c("Screening", "Diagnostic testing"))) |>
    filter(year %in% c(2025:2030)) |>
    group_by(lhs_num, test_pop, test_strat, test_num, year) |>
    summarise(n_I = sum(n_I), NS = sum(NS), I = sum(I)) |>
    ungroup() |>
    mutate(incidence_intervention = n_I/(NS - I)) |>
    group_by(lhs_num, test_pop, test_strat, year) |>
    # Calculate percentage change relative to baseline (test_num = 0)
    mutate(incidence_baseline = incidence_intervention[test_num == 0],
           prop_incidence_reduction = (incidence_baseline - incidence_intervention) / incidence_baseline) |>
    select(lhs_num, test_strat, test_pop, test_num, incidence_baseline, incidence_intervention, prop_incidence_reduction) |>
    ungroup() 
)

fwrite(df_year, "./intervention/impact_intervention_constrained_year.csv")

# SENSITIVITY ANALYSIS ----

lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

df_intervention_strata <- data.table::fread("./posterior/post_intervention_constrained_strata_sensitivity.csv") |>
  filter(lhs_num %in% lhsnum)

df_intervention <- df_intervention_strata[year %in% 2024:2030,
                                          lapply(.SD, sum),
                                          by = .(lhs_num, year, test_pop, test_strat, test_num),
                                          .SDcols = !c("sex", "risk", "age", "gest", "sex2", "risk2", "age2", "gest2")]

df_intervention_total <- df_intervention |>
  mutate(test_strat = factor(test_strat, levels = c("h","y"),
                             labels = c("Screening", "Diagnostic testing"))) |>
  mutate(eligible = case_when(test_strat == "Screening" ~ n_elig_h,
                              test_strat == "Diagnostic testing" ~ n_elig_y + n_elig_o),
         test_used = case_when(test_strat == "Screening" ~ n_test_h,
                               test_strat == "Diagnostic testing" ~ n_test_y + n_test_o))

df_intervention_total |> summarise(n = n_distinct(lhs_num))

# Calculate outcomes:
# cases averted and prop cases averted relative to baseline, 2025-2030
# incidence reduction and prop incidence reduction relative to baseline, 2030
df <- left_join(
  df_intervention_total |> 
    filter(year %in% c(2025:2030)) |>
    group_by(lhs_num, test_pop, test_strat, test_num) |>
    summarise(cases = sum(n_I),
              eligible = sum(eligible),
              test_used = sum(test_used)) |>
    ungroup() |>
    group_by(lhs_num, test_pop, test_strat) |>
    # Calculate change relative to baseline (test_num = 0)
    mutate(cases_baseline = cases[test_num == 0],
           cases_averted = cases_baseline - cases,
           prop_cases_averted = cases_averted / cases_baseline,
           test_avail = test_num*6) |>
    ungroup() |>
    filter(test_num != 0) |>
    rename(cases_intervention = cases) |>
    relocate(lhs_num, test_strat, test_pop, test_num, eligible, test_avail, test_used, cases_baseline, cases_intervention, cases_averted, prop_cases_averted),
  
  df_intervention_total |> 
    filter(year == 2030) |>
    mutate(incidence_intervention = n_I/(NS - I)) |>
    group_by(lhs_num, test_pop, test_strat) |>
    # Calculate percentage change relative to baseline (test_num = 0)
    mutate(incidence_baseline = incidence_intervention[test_num == 0],
           prop_incidence_reduction = (incidence_baseline - incidence_intervention) / incidence_baseline) |>
    select(lhs_num, test_strat, test_pop, test_num, incidence_baseline, incidence_intervention, prop_incidence_reduction) |>
    ungroup() |>
    filter(test_num != 0)
)

write.csv(df, "./intervention/impact_intervention_constrained_total_sensitivity.csv", row.names = FALSE)

# FOR MAX NUMBER POCT ----

post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)

lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

df_intervention_max <- data.table::fread("./posterior/post_intervention_unrestricted_strata.csv") |>
  filter(lhs_num %in% lhsnum)

df_intervention_max |> group_by() |> summarise(n=n_distinct(lhs_num))
df_intervention_max |> count(test_strat, test_pop, test_num)

# by sex
df_intervention_max_sex <- df_intervention_max[year %in% 2024:2030,
                                              lapply(.SD, sum),
                                              by = .(lhs_num, year, test_pop, test_strat, test_num, sex),
                                              .SDcols = !c("risk", "age", "gest", "sex2", "risk2", "age2", "gest2")]

df_intervention_max_sex <- df_intervention_max_sex |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
         test_strat = factor(test_strat, levels = c("h","y"),
                             labels = c("Screening", "Diagnostic testing"))) |>
  mutate(eligible = case_when(test_strat == "Screening" ~ n_elig_h,
                              test_strat == "Diagnostic testing" ~ n_elig_y + n_elig_o),
         test_used = case_when(test_strat == "Screening" ~ n_test_h,
                               test_strat == "Diagnostic testing" ~ n_test_y + n_test_o))

fwrite(df_intervention_max_sex, "./posterior/post_intervention_unrestricted_sex.csv")

# overall
df_intervention_max_total <- df_intervention_max[year %in% 2024:2030,
                                          lapply(.SD, sum),
                                          by = .(lhs_num, year, test_pop, test_strat, test_num),
                                          .SDcols = !c("sex", "risk", "age", "gest", "sex2", "risk2", "age2", "gest2")]

df_intervention_max_total <- df_intervention_max_total |>
  mutate(test_strat = factor(test_strat, levels = c("h","y"),
                             labels = c("Screening", "Diagnostic testing"))) |>
  mutate(eligible = case_when(test_strat == "Screening" ~ n_elig_h,
                              test_strat == "Diagnostic testing" ~ n_elig_y + n_elig_o),
         test_used = case_when(test_strat == "Screening" ~ n_test_h,
                               test_strat == "Diagnostic testing" ~ n_test_y + n_test_o))

fwrite(df_intervention_max_total, "./posterior/post_intervention_unrestricted_total.csv")

# IMPACT ----

## Overall----

# Import data if has already been defined and saved above
df_intervention_max_total <- read.csv("./posterior/post_intervention_unrestricted_total.csv") 

df_intervention_max_total |> summarise(n = n_distinct(lhs_num))

# Calculate outcomes:
# cases averted and prop cases averted relative to baseline, 2025-2030
# incidence reduction and prop incidence reduction relative to baseline, 2030
df <- left_join(
  df_intervention_max_total |> 
    filter(year %in% c(2025:2030)) |>
    group_by(lhs_num, test_pop, test_strat, test_num) |>
    summarise(cases = sum(n_I),
              eligible = sum(eligible),
              test_used = sum(test_used)) |>
    ungroup() |>
    group_by(lhs_num, test_pop, test_strat) |>
    # Calculate change relative to baseline (test_num = 0)
    mutate(cases_baseline = cases[test_num == 0],
           cases_averted = cases_baseline - cases,
           prop_cases_averted = cases_averted / cases_baseline) |>
    ungroup() |>
    filter(test_num != 0) |>
    rename(cases_intervention = cases) |>
    relocate(lhs_num, test_strat, test_pop, test_num, eligible, test_used, cases_baseline, cases_intervention, cases_averted, prop_cases_averted),
  
  df_intervention_max_total |> 
    filter(year == 2030) |>
    mutate(incidence_intervention = n_I/(NS - I)) |>
    group_by(lhs_num, test_pop, test_strat) |>
    # Calculate percentage change relative to baseline (test_num = 0)
    mutate(incidence_baseline = incidence_intervention[test_num == 0],
           prop_incidence_reduction = (incidence_baseline - incidence_intervention) / incidence_baseline) |>
    select(lhs_num, test_strat, test_pop, test_num, incidence_baseline, incidence_intervention, prop_incidence_reduction) |>
    ungroup() |>
    filter(test_num != 0)
)

write.csv(df, "./intervention/impact_intervention_unrestricted_total.csv", row.names = FALSE)



df_intervention_max_total |> 
  filter(year %in% c(2024:2030)) |>
  group_by(lhs_num, test_pop, test_strat, test_num, year) |>
  summarise(cases = sum(n_I),
            eligible = sum(eligible),
            test_used = sum(test_used)) |>
  ungroup() |>
  ggplot() +
  geom_boxplot(aes(x = as.factor(year), y = cases, colour = as.factor(test_num))) +
  facet_grid(test_strat~test_pop)
  
  group_by(lhs_num, test_pop, test_strat) |>
  # Calculate change relative to baseline (test_num = 0)
  mutate(cases_baseline = cases[test_num == 0],
         cases_averted = cases_baseline - cases,
         prop_cases_averted = cases_averted / cases_baseline) |>
  ungroup() 
