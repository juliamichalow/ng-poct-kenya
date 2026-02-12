# Calculations for baseline scenario

lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

intervention_strata <- data.table::fread("./posterior/post_intervention_constrained_strata.csv") |>
  filter(lhs_num %in% lhsnum)

inf_strata <- intervention_strata |>
  filter(test_num == 0, test_pop == "fsw", test_strat == "h",
         year %in% c(2025:2030)) |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
         risk = factor(risk, levels = c(1,2,3), labels = c("Low", "Medium", "High")),
         age = factor(age, levels = c(1,2), labels = c("Young", "Old")),
         gest = factor(gest, levels = c(1,2), labels = c("Non-pregnant", "Pregnant"))) |>
  filter(!(sex == "Male" & gest == "Pregnant"))

rm(intervention_strata)

# Prev ----

inf_strata |>
  filter(year == 2025) |> 
  group_by(lhs_num, sex) |>
  summarise(value = sum(I)/sum(NS)) |>
  group_by(sex) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))

calc_prev <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(I)/sum(NS)) |>
    group_by() |>
    summarise(median = median(value), 
              lwr = quantile(value, 0.025), 
              upr = quantile(value, 0.975))
}

table_prev <- rbind(
  calc_prev(inf_strata, year == 2025, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_prev(inf_strata, year == 2025, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "PREG"),
  calc_prev(inf_strata, year == 2025, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_prev(inf_strata, year == 2025, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "MEN"),
  calc_prev(inf_strata, year == 2025, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("PREG","AGYW","FSW","MEN", "CFSW")))

# Incidence ----

# Incident cases total
inf_strata |>
  filter(year == 2025) |> 
  group_by(lhs_num) |>
  summarise(value = sum(n_I)) |>
  group_by() |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))

# Prop incident cases high risk
inf_strata |>
  filter(year == 2025) |>
  group_by(lhs_num, risk) |>
  summarise(value = sum(n_I)) |>
  pivot_wider(names_from = risk, values_from = value) |>
  mutate(value = High / (Low + Medium + High)) |>
  group_by() |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 


# Prop incident cases high risk by sex
inf_strata |>
  filter(year == 2025) |>
  group_by(lhs_num, sex, risk) |>
  summarise(value = sum(n_I)) |>
  pivot_wider(names_from = risk, values_from = value) |>
  mutate(value = High / (Low + Medium + High)) |>
  group_by(sex) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 


# Prop incident cases male
inf_strata |>
  filter(year == 2025) |>
  group_by(lhs_num, sex) |>
  summarise(value = sum(n_I)) |>
  pivot_wider(names_from = sex, values_from = value) |>
  mutate(value = Male / (Male + Female)) |>
  group_by() |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))


# Incident cases per 100 population
inf_strata |>
  filter(year == 2025) |>
  group_by(lhs_num, sex) |>
  summarise(value = sum(n_I)/sum(NS)) |>
  group_by(sex) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 


inf_strata |>
  filter(year == 2025, risk == "High") |>
  group_by(lhs_num, sex) |>
  summarise(value = sum(n_I)/sum(NS)) |>
  group_by(sex) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 


calc_inc <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(n_I)/(sum(NS) - sum(I))) |>
    group_by() |>
    summarise(median = median(value), 
              lwr = quantile(value, 0.025), 
              upr = quantile(value, 0.975))
}

table_inc <- rbind(
  calc_inc(inf_strata, year == 2025, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_inc(inf_strata, year == 2025, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "PREG"),
  calc_inc(inf_strata, year == 2025, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_inc(inf_strata, year == 2025, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "MEN"),
  calc_inc(inf_strata, year == 2025, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("PREG","AGYW","FSW","MEN", "CFSW")))

calc_case <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(n_I)) |>
    group_by() |>
    summarise(median = median(value), 
              lwr = quantile(value, 0.025), 
              upr = quantile(value, 0.975))
}

table_case <- rbind(
  calc_case(inf_strata, year == 2025, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_case(inf_strata, year == 2025, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "PREG"),
  calc_case(inf_strata, year == 2025, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_case(inf_strata, year == 2025, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "MEN"),
  calc_case(inf_strata, year == 2025, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("PREG","AGYW","FSW","MEN", "CFSW")))

# Symptoms and treatment ----

inf_strata |>
  filter(year == 2025) |>
  group_by(lhs_num, sex) |>
  summarise(prop_symp = (sum(n_Y) + sum(n_Z))/sum(n_I),
            prop_tx_symp = sum(n_Y)/(sum(n_Y) + sum(n_Z)),
            ap = sum(n_Y)/sum(n_symp),
            n_symp = sum(n_symp)) |>
  pivot_longer(cols = c(prop_symp, prop_tx_symp, ap, n_symp)) |>
  group_by(sex,name) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))
  
  
# QALY baseline ----

qaly_strata <- read.csv("../data/impact_qaly_tot_strata_2025.csv") |>
  filter(test_pop == "agyw_sa",
         test_strat == "Screening") |>
  select(!c(test_pop, test_strat, test_num, test_used))

calc_qaly <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(qaly_baseline)) |>
    group_by() |>
    summarise(median = median(value), 
              lwr = quantile(value, 0.025), 
              upr = quantile(value, 0.975))
}

table_qaly <- rbind(
  calc_qaly(qaly_strata, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_qaly(qaly_strata, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "PREG"),
  calc_qaly(qaly_strata, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_qaly(qaly_strata, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "MEN"),
  calc_qaly(qaly_strata, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("PREG","AGYW","FSW","MEN", "CFSW")))

# QALYs lost total
qaly_strata |>
  group_by(lhs_num) |>
  summarise(value = sum(qaly_baseline)) |>
  group_by() |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))


# Prop QALYs lost high risk
qaly_strata |>
  group_by(lhs_num, risk) |>
  summarise(value = sum(qaly_baseline)) |>
  pivot_wider(names_from = risk, values_from = value) |>
  mutate(value = High / (Low + Medium + High)) |>
  group_by() |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))

# Prop QALYs lost pregnant
qaly_strata |>
  group_by(lhs_num, sex, gest) |>
  summarise(value = sum(qaly_baseline)) |>
  pivot_wider(names_from = c(sex,gest), values_from = value) |>
  mutate(tot = `Female_Non-pregnant` + Female_Pregnant + `Male_Non-pregnant`) |>
  mutate(prop_m = `Male_Non-pregnant` / tot,
         prop_f_non = `Female_Non-pregnant` / tot,
         prop_f_preg = Female_Pregnant / tot) |>
  select(lhs_num, starts_with("prop")) |>
  pivot_longer(cols = c("prop_m", "prop_f_non", "prop_f_preg")) |>
  group_by(name) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))

# Pop ----

# Prop pop high risk
inf_strata |>
  filter(year == 2025) |>
  group_by(lhs_num, risk) |>
  summarise(value = sum(N)) |>
  pivot_wider(names_from = risk, values_from = value) |>
  mutate(value = High / (Low + Medium + High)) |>
  group_by() |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 


inf_strata |>
  filter(year == 2025) |>
  group_by(lhs_num, sex, risk) |>
  summarise(value = sum(NS)) |>
  pivot_wider(names_from = risk, values_from = value) |>
  mutate(value = High / (Low + Medium + High)) |>
  group_by(sex) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 

# Percentage high risk of total pop
inf_strata |>
  filter(year == 2025, risk == "High") |>
  group_by(lhs_num, sex, risk) |>
  summarise(pop = sum(NS)) |>
  left_join(inf_strata |>
              filter(year == 2025) |>
              group_by(lhs_num) |>
              summarise(tot = sum(NS))) |>
  mutate(value = pop/tot) |>
  group_by(sex) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 
