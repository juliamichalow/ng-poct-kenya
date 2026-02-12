# MINIMUM SCENARIO ----

# INPUT DATA ----

# Baseline 2025

lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

intervention_strata <- data.table::fread("./posterior/post_intervention_constrained_strata.csv") |>
  filter(lhs_num %in% lhsnum)

inf_strata <- intervention_strata |>
  filter(test_num == 0, test_pop == "fsw", test_strat == "h",
         year == 2025) |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
         risk = factor(risk, levels = c(1,2,3), labels = c("Low", "Medium", "High")),
         age = factor(age, levels = c(1,2), labels = c("Young", "Old")),
         gest = factor(gest, levels = c(1,2), labels = c("Non-pregnant", "Pregnant")))

rm(intervention_strata)

qaly_strata <- read.csv("./qalys/impact_qaly_constrained_tot_strata_2025.csv") |>
  filter(test_pop == "agyw_sa",
         test_strat == "Screening") |>
  select(!c(test_pop, test_strat, test_num, test_used))

# Results summary

impact_cases_year <- read.csv("./posterior/post_intervention_constrained_total.csv") |>
  mutate(target_pop = factor(test_pop, levels = c("fsw","male_h","agyw_sa","preg","male_all"),
                             labels = c("FSW","CFSW","AGYW","PREG","MEN")),
         test_strat = factor(test_strat, levels = c("Screening", "Diagnostic testing"))) |>
  filter(year %in% c(2025:2030))

impact_case <- read.csv("./intervention/impact_intervention_constrained_total.csv") |>
  filter(test_num == 25000) |>
  mutate(target_pop = factor(test_pop, levels = c("fsw","male_h","agyw_sa","preg","male_all"),
                             labels = c("FSW","CFSW","AGYW","PREG","MEN")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))

impact_qaly <- read.csv("./qalys/impact_qaly_constrained_tot.csv") |>
  filter(test_num == 25000) |>
  mutate(target_pop = factor(test_pop, levels = c("fsw","male_h","agyw_sa","preg","male_all"),
                             labels = c("FSW","CFSW","AGYW","PREG","MEN")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))

# PREP DATA ----

calc_base <- function(inf, qaly, ...) {
  
  rbind(inf |>
          filter(...) |>
          group_by(lhs_num) |>
          summarise(N = sum(N),
                    NS = sum(NS),
                    inc = sum(n_I)/(sum(NS) - sum(I)),
                    prev = sum(I)/sum(NS)) |>
          pivot_longer(cols = c(N, NS, inc, prev), names_to = "var") |>
          group_by(var) |>
          summarise(median = median(value, na.rm=TRUE), 
                    lwr = quantile(value, 0.025, na.rm=TRUE), 
                    upr = quantile(value, 0.975, na.rm=TRUE)),
        qaly |>
          filter(...) |>
          group_by(lhs_num) |>
          summarise(qaly = sum(qaly_baseline)) |>
          left_join(qaly |> 
                      group_by(lhs_num) |> 
                      summarise(tot = sum(qaly_baseline))) |>
          mutate(qaly_oftot = qaly/tot) |>
          pivot_longer(cols = c(qaly,qaly_oftot), names_to = "var") |>
          group_by(var) |>
          summarise(median = median(value, na.rm=TRUE), 
                    lwr = quantile(value, 0.025, na.rm=TRUE), 
                    upr = quantile(value, 0.975, na.rm=TRUE)),
        
        inf |>
          filter(...) |>
          group_by(lhs_num) |>
          summarise(case = sum(n_I)) |>
          left_join(inf |> 
                      group_by(lhs_num) |> 
                      summarise(tot = sum(n_I))) |>
          mutate(case_oftot = case/tot) |>
          pivot_longer(cols = case_oftot, names_to = "var") |>
          group_by(var) |>
          summarise(median = median(value, na.rm=TRUE), 
                    lwr = quantile(value, 0.025, na.rm=TRUE), 
                    upr = quantile(value, 0.975, na.rm=TRUE)))
  }

dat_base <- rbind(
  calc_base(inf_strata, qaly_strata, 
       sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_base(inf_strata, qaly_strata, 
       sex == "Female", gest == "Pregnant") |> mutate(target_pop = "PREG"),
  calc_base(inf_strata, qaly_strata, 
       sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_base(inf_strata, qaly_strata, 
       sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "MEN"),
  calc_base(inf_strata, qaly_strata, 
       sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) 

dat_int <- 
  impact_cases_year |>
      filter(year == 2025) |>
      filter(test_num == 0) |>
      group_by(test_strat, target_pop) |>
      summarise(median = median(eligible), 
                lwr = quantile(eligible, 0.025), 
                upr = quantile(eligible, 0.975), .groups = "drop") |>
      mutate(var= case_when(test_strat == "Screening" ~ "eligible_screen",
                            test_strat == "Diagnostic testing" ~ "eligible_diag")) |>
      select(!test_strat) |>
  full_join(
    impact_cases_year |>
      filter(test_num == 25000) |>
      mutate(prop_test = test_used/eligible) |>
      group_by(test_strat, target_pop) |>
      summarise(median = median(prop_test),
                lwr = quantile(prop_test, 0.025), 
                upr = quantile(prop_test, 0.975), .groups = "drop") |>
      mutate(var= case_when(test_strat == "Screening" ~ "prop_screen",
                            test_strat == "Diagnostic testing" ~ "prop_diag")) |>
      select(!test_strat)) |>
  full_join(
    impact_case |> 
      mutate(pertest_cases = cases_averted/test_used, prop_cases = prop_cases_averted) |>
      select(lhs_num, test_strat, target_pop, prop_cases, pertest_cases) |>
      pivot_longer(cols = c(prop_cases, pertest_cases)) |>
      group_by(test_strat, target_pop, name) |>
      summarise(median = median(value, na.rm=TRUE), 
                lwr = quantile(value, 0.025, na.rm=TRUE), 
                upr = quantile(value, 0.975, na.rm=TRUE),.groups = "drop") |>
      mutate(var = case_when(test_strat == "Screening" & name == "prop_cases" ~ "caseavert_screen",
                             test_strat == "Diagnostic testing"  & name == "prop_cases"  ~ "caseavert_diag",
                             test_strat == "Screening" & name == "pertest_cases" ~ "casepertest_screen",
                             test_strat == "Diagnostic testing"  & name == "pertest_cases"  ~ "casepertest_diag")) |> 
      select(!c(test_strat,name))) |>
  full_join(
    impact_qaly |> 
      mutate(prop_qaly = qaly_gained/qaly_baseline, pertest_qaly = qaly_gained/test_used) |>
      select(lhs_num, test_strat, target_pop, prop_qaly, pertest_qaly) |>
      pivot_longer(cols = c(prop_qaly, pertest_qaly)) |>
      group_by(test_strat, target_pop, name) |>
      summarise(median = median(value, na.rm=TRUE), 
                lwr = quantile(value, 0.025, na.rm=TRUE), 
                upr = quantile(value, 0.975, na.rm=TRUE),.groups = "drop") |>
      mutate(var = case_when(test_strat == "Screening" & name == "prop_qaly" ~ "qalygain_screen",
                             test_strat == "Diagnostic testing"  & name == "prop_qaly"  ~ "qalygain_diag",
                             test_strat == "Screening" & name == "pertest_qaly" ~ "qalypertest_screen",
                             test_strat == "Diagnostic testing"  & name == "pertest_qaly"  ~ "qalypertest_diag")) |>
      select(!c(test_strat,name))) 

dat <- full_join(dat_base, dat_int) 

tab <- dat |>
  pivot_longer(cols = c(median, lwr, upr)) |>
  pivot_wider(names_from = var, values_from = value) |>
  mutate(qaly_perpop = qaly/NS) |>
  select(!c(qaly, NS)) |>
  mutate(N = case_when(N>10*1000000 ~ sprintf("%.1f", N/1000000), TRUE ~ sprintf("%.2f", N/1000000)),
         prev = case_when(prev>0.1 ~ sprintf("%.1f", prev*100), TRUE ~ sprintf("%.2f", prev*100)),
         inc = case_when(inc>0.1 ~ sprintf("%.1f", inc*100), TRUE ~ sprintf("%.2f", inc*100)),
         qaly_perpop = case_when(qaly_perpop>0.1 ~ sprintf("%.0f", qaly_perpop*1000),
                                 qaly_perpop>0.01 ~ sprintf("%.1f", qaly_perpop*1000), TRUE ~ sprintf("%.2f", qaly_perpop*1000)),
         qaly_oftot = case_when(qaly_oftot>0.1 ~ sprintf("%.1f", qaly_oftot*100), TRUE ~ sprintf("%.2f", qaly_oftot*100)),
         case_oftot = case_when(case_oftot>0.1 ~ sprintf("%.1f", case_oftot*100), TRUE ~ sprintf("%.2f", case_oftot*100)),
         prop_screen = case_when(prop_screen < 0.1 ~ sprintf("%.2f", prop_screen*100), TRUE ~ sprintf("%.1f", prop_screen*100)),
         prop_diag = case_when(prop_diag < 0.1 ~ sprintf("%.2f", prop_diag*100), prop_diag >= 1 ~ sprintf("%.0f", prop_diag*100),
                               TRUE ~ sprintf("%.1f", prop_diag*100)),
         caseavert_diag = sprintf("%.2f", caseavert_diag*100),
         caseavert_screen = sprintf("%.2f", caseavert_screen*100),
         qalygain_diag = sprintf("%.2f", qalygain_diag*100),
         qalygain_screen = sprintf("%.2f", qalygain_screen*100),
         casepertest_diag = case_when(casepertest_diag>1 ~ sprintf("%.0f", casepertest_diag*100), 
                                      casepertest_diag>0.1 ~ sprintf("%.1f", casepertest_diag*100), 
                                      TRUE ~ sprintf("%.2f", casepertest_diag*100)),
         casepertest_screen = case_when(casepertest_screen>1 ~ sprintf("%.0f", casepertest_screen*100), 
                                        casepertest_screen>0.1 ~ sprintf("%.1f", casepertest_screen*100), 
                                        TRUE ~ sprintf("%.2f", casepertest_screen*100)),
         qalypertest_diag = case_when(qalypertest_diag>1 ~ sprintf("%.0f", qalypertest_diag*100), 
                                      qalypertest_diag>0.1 ~ sprintf("%.1f", qalypertest_diag*100), 
                                      TRUE ~ sprintf("%.2f", qalypertest_diag*100)),
         qalypertest_screen = case_when(qalypertest_screen>1 ~ sprintf("%.0f", qalypertest_screen*100), 
                                        qalypertest_screen>0.1 ~ sprintf("%.1f", qalypertest_screen*100), 
                                        TRUE ~ sprintf("%.2f", qalypertest_screen*100)),
         eligible_screen = case_when(eligible_screen <= 10*1e6 ~ sprintf("%.2f",eligible_screen/1000000), TRUE ~ sprintf("%.2f",eligible_screen/1000000)),
         eligible_diag = case_when(eligible_diag <= 10*1e6 ~ sprintf("%.2f",eligible_diag/1000000), TRUE ~ sprintf("%.2f",eligible_diag/1000000))) |>
  pivot_longer(cols = !c(target_pop, name), names_to = "var") |>
  pivot_wider(names_from = name, values_from = value) |>
  mutate(est = paste0(median, " (", lwr, " - ", upr, ")")) |>
  select(!c(median, lwr, upr)) |>
  pivot_wider(names_from = target_pop, values_from = est) |>
  mutate(var= factor(var, levels = c("N","prev", "inc", "qaly_perpop","case_oftot","qaly_oftot", 
                                     "eligible_diag", "prop_diag", "caseavert_diag", "qalygain_diag", "casepertest_diag", "qalypertest_diag",
                                     "eligible_screen", "prop_screen", "caseavert_screen", "qalygain_screen", "casepertest_screen", "qalypertest_screen"),
                     labels = c("Population size (millions)","Prevalence (%)", "Incidence (per 100)", "QALYs lost (per 1000)",
                                "Proportion total incident cases (%)", "Proportion total QALYs lost (%)", 
                                "Eligible consultations d (millions)", "Eligible consultations tested d (%)", "Cases averted relative to baseline d (%)", "QALYs gained relative to baseline d (%)", "Cases averted per 100 tests d", "QALYs gained per 100 tests d",
                                "Eligible consultations s (millions)", "Eligible consultations tested s (%)", "Cases averted relative to baseline s (%)", "QALYs gained relative to baseline s (%)", "Cases averted per 100 tests s", "QALYs gained per 100 tests s"))) |>
  arrange(var) |>
  mutate(latex = paste0(CFSW," & ", FSW," & ",MEN," & ",PREG," & ",AGYW, "\\")) 

write.csv(tab, ".././tables/results_baseline.csv")


# MAXIMUM SCENARIO ----

# INPUT DATA ----

# Results summary

impact_cases_year <- read.csv("./posterior/post_intervention_unrestricted_total.csv") |>
  mutate(target_pop = factor(test_pop, levels = c("fsw","male_h","agyw_sa","preg","male_all"),
                             labels = c("FSW","CFSW","AGYW","PREG","MEN")),
         test_strat = factor(test_strat, levels = c("Screening", "Diagnostic testing"))) |>
  filter(year %in% c(2025:2030))

impact_case <- read.csv("./intervention/impact_intervention_unrestricted_total.csv") |>
  filter(test_num== 10) |>
  mutate(target_pop = factor(test_pop, levels = c("fsw","male_h","agyw_sa","preg","male_all"),
                             labels = c("FSW","CFSW","AGYW","PREG","MEN")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))

impact_qaly <- read.csv("./qalys/impact_qaly_unrestricted_tot.csv") |>
  filter(test_num == 10) |>
  mutate(target_pop = factor(test_pop, levels = c("fsw","male_h","agyw_sa","preg","male_all"),
                             labels = c("FSW","CFSW","AGYW","PREG","MEN")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))

# PREP DATA ----

dat_int <- 
  impact_cases_year |>
  filter(test_num == 0) |>
  group_by(lhs_num, test_strat, target_pop) |>
  summarise(eligible = sum(eligible)/6) |>
  ungroup() |>
  group_by(test_strat, target_pop) |>
  summarise(median = median(eligible), 
            lwr = quantile(eligible, 0.025), 
            upr = quantile(eligible, 0.975), .groups = "drop") |>
  mutate(var= case_when(test_strat == "Screening" ~ "eligible_screen",
                        test_strat == "Diagnostic testing" ~ "eligible_diag")) |>
  select(!test_strat) |>
  full_join(
    impact_cases_year |>
      filter(test_num == 10) |>
      group_by(lhs_num, test_strat, target_pop) |>
      mutate(prop_test = sum(test_used)/sum(eligible)) |>
      ungroup() |>
      group_by(test_strat, target_pop) |>
      summarise(median = median(prop_test),
                lwr = quantile(prop_test, 0.025), 
                upr = quantile(prop_test, 0.975), .groups = "drop") |>
      mutate(var= case_when(test_strat == "Screening" ~ "prop_screen",
                            test_strat == "Diagnostic testing" ~ "prop_diag")) |>
      select(!test_strat)) |>
  full_join(
    impact_case |> 
      mutate(pertest_cases = cases_averted/test_used, prop_cases = prop_cases_averted) |>
      select(lhs_num, test_strat, target_pop, prop_cases, pertest_cases) |>
      pivot_longer(cols = c(prop_cases, pertest_cases)) |>
      group_by(test_strat, target_pop, name) |>
      summarise(median = median(value, na.rm=TRUE), 
                lwr = quantile(value, 0.025, na.rm=TRUE), 
                upr = quantile(value, 0.975, na.rm=TRUE),.groups = "drop") |>
      mutate(var = case_when(test_strat == "Screening" & name == "prop_cases" ~ "caseavert_screen",
                             test_strat == "Diagnostic testing"  & name == "prop_cases"  ~ "caseavert_diag",
                             test_strat == "Screening" & name == "pertest_cases" ~ "casepertest_screen",
                             test_strat == "Diagnostic testing"  & name == "pertest_cases"  ~ "casepertest_diag")) |> 
      select(!c(test_strat,name))) |>
  full_join(
    impact_qaly |> 
      mutate(prop_qaly = qaly_gained/qaly_baseline, pertest_qaly = qaly_gained/test_used) |>
      select(lhs_num, test_strat, target_pop, prop_qaly, pertest_qaly) |>
      pivot_longer(cols = c(prop_qaly, pertest_qaly)) |>
      group_by(test_strat, target_pop, name) |>
      summarise(median = median(value, na.rm=TRUE), 
                lwr = quantile(value, 0.025, na.rm=TRUE), 
                upr = quantile(value, 0.975, na.rm=TRUE),.groups = "drop") |>
      mutate(var = case_when(test_strat == "Screening" & name == "prop_qaly" ~ "qalygain_screen",
                             test_strat == "Diagnostic testing"  & name == "prop_qaly"  ~ "qalygain_diag",
                             test_strat == "Screening" & name == "pertest_qaly" ~ "qalypertest_screen",
                             test_strat == "Diagnostic testing"  & name == "pertest_qaly"  ~ "qalypertest_diag")) |>
      select(!c(test_strat,name))) 

tab <- dat_int |>
  pivot_longer(cols = c(median, lwr, upr)) |>
  pivot_wider(names_from = var, values_from = value) |>
  mutate(caseavert_diag = sprintf("%.2f", caseavert_diag*100),
         caseavert_screen = sprintf("%.2f", caseavert_screen*100),
         qalygain_diag = sprintf("%.2f", qalygain_diag*100),
         qalygain_screen = sprintf("%.2f", qalygain_screen*100),
         casepertest_diag = case_when(casepertest_diag>1 ~ sprintf("%.0f", casepertest_diag*100), 
                                      casepertest_diag>0.1 ~ sprintf("%.1f", casepertest_diag*100), 
                                      TRUE ~ sprintf("%.2f", casepertest_diag*100)),
         casepertest_screen = case_when(casepertest_screen>1 ~ sprintf("%.0f", casepertest_screen*100), 
                                        casepertest_screen>0.1 ~ sprintf("%.1f", casepertest_screen*100), 
                                        TRUE ~ sprintf("%.2f", casepertest_screen*100)),
         qalypertest_diag = case_when(qalypertest_diag>1 ~ sprintf("%.0f", qalypertest_diag*100), 
                                      qalypertest_diag>0.1 ~ sprintf("%.1f", qalypertest_diag*100), 
                                      TRUE ~ sprintf("%.2f", qalypertest_diag*100)),
         qalypertest_screen = case_when(qalypertest_screen>1 ~ sprintf("%.0f", qalypertest_screen*100), 
                                        qalypertest_screen>0.1 ~ sprintf("%.1f", qalypertest_screen*100), 
                                        TRUE ~ sprintf("%.2f", qalypertest_screen*100)),
         eligible_screen = sprintf("%.2f",eligible_screen/100000),
         eligible_diag = sprintf("%.2f",eligible_diag/100000)) |>
  select(!c(prop_diag, prop_screen)) |>
  pivot_longer(cols = !c(target_pop, name), names_to = "var") |>
  pivot_wider(names_from = name, values_from = value) |>
  mutate(est = paste0(median, " (", lwr, " - ", upr, ")")) |>
  select(!c(median, lwr, upr)) |>
  pivot_wider(names_from = target_pop, values_from = est) |>
  mutate(var= factor(var, levels = c("eligible_diag", "prop_diag", "caseavert_diag", "qalygain_diag", "casepertest_diag", "qalypertest_diag",
                                     "eligible_screen", "prop_screen", "caseavert_screen", "qalygain_screen", "casepertest_screen", "qalypertest_screen"),
                     labels = c("Eligible consultations d (100000s)", "Eligible consultations tested d (%)", "Cases averted relative to baseline d (%)", "QALYs gained relative to baseline d (%)", "Cases averted per 100 tests d", "QALYs gained per 100 tests d",
                                "Eligible consultations s (100000s)", "Eligible consultations tested s (%)", "Cases averted relative to baseline s (%)", "QALYs gained relative to baseline s (%)", "Cases averted per 100 tests s", "QALYs gained per 100 tests s"))) |>
  arrange(var) |>
  mutate(latex = paste0(CFSW," & ", FSW," & ",MEN," & ",PREG," & ",AGYW, "\\")) 

tab |> write.csv(".././tables/results_unrestricted.csv")
  