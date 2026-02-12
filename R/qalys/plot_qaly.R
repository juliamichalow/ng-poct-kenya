library(tidyverse)
library(patchwork)
library(ggh4x)

# INPUT DATA
# Stratified by sequelae
qaly_seq_sexgest <- read.csv("./qalys/impact_qaly_constrained_seq_sexgest.csv") 
qaly_seq_strata_2025 <- read.csv("./qalys/impact_qaly_constrained_seq_strata_2025.csv") 

# Aggregated across sequelae
qaly_tot_sexgest <- read.csv("./qalys/impact_qaly_constrained_tot_sexgest.csv") 
qaly_tot_strata <- read.csv("./qalys/impact_qaly_constrained_tot_strata.csv") 
qaly_tot_strata_2025 <- read.csv("./qalys/impact_qaly_constrained_tot_strata_2025.csv") 

mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=-10,r=5,b=0,l=5), "cm")),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

# ANALYSE ----

## Baseline results ----

## QALYs lost per 1000 incident infections
qaly_tot_strata_2025 |>
  filter(test_pop == "agyw_sa", test_strat == "Screening", test_num == 0) |>
  group_by(lhs_num, sex, gest) |>
  summarise(case = sum(case), qaly = sum(qaly_baseline)) |>
  mutate(qaly_per_case = qaly / case * 1000) |>
  ungroup() |>
  group_by(sex, gest) |>
  summarise(median = median(qaly_per_case),
            lwr = quantile(qaly_per_case, 0.025),
            upr = quantile(qaly_per_case, 0.975))

## QALYS lost per person-years
qaly_tot_strata_2025 |>
  filter(test_pop == "agyw_sa", test_strat == "Screening", test_num == 0) |>
  group_by(lhs_num, sex, gest) |>
  summarise(py = sum(annual_py), qaly = sum(qaly_baseline)) |>
  mutate(qaly_per_py = qaly / py * 1000) |>
  ungroup() |>
  group_by(sex, gest) |>
  summarise(median = median(qaly_per_py),
            lwr = quantile(qaly_per_py, 0.025),
            upr = quantile(qaly_per_py, 0.975))

## Population-level discounted QALYs lost 
qaly_tot_strata_2025 |>
  filter(test_pop == "agyw_sa", test_strat == "Screening", test_num == 0) |>
  group_by(lhs_num, sex, gest) |>
  summarise(qaly = sum(qaly_baseline)) |>
  ungroup() |>
  group_by(sex, gest) |>
  summarise(median = median(qaly),
            lwr = quantile(qaly, 0.025),
            upr = quantile(qaly, 0.975))

## Composition of QALYs lost
qaly_tot_strata_2025 |>
  filter(test_strat == "Screening", test_pop == "agyw_sa", test_num == 0) |>
  group_by(lhs_num, sex, gest) |>
  summarise(qaly_baseline = sum(qaly_baseline), .groups = "drop") |>
  group_by(lhs_num) |>
  mutate(prop = qaly_baseline/sum(qaly_baseline)) |>
  ungroup() |>
  group_by(sex, gest) |>
  summarise(median = median(prop),
            lwr = quantile(prop, 0.025),
            upr = quantile(prop, 0.975))

# PLOTS ----

# Composition of QALYs lost
df_plot <- 
  qaly_seq_strata_2025 |>
  filter(test_strat == "Screening", test_pop == "agyw_sa") |>
  pivot_longer(cols = c(starts_with("qaly")), names_to = "var") |>
  mutate(var = str_replace(var, "qaly_", ""),
         pop = case_when(sex == "Male" ~ "Male", 
                         sex == "Female" & gest == "Non-pregnant" ~ "Female, Non-pregnant",
                         sex == "Female" & gest == "Pregnant" & var %in% c("n_still", "n_lbw", "n_pneu", "n_conj") ~ "Female, Pregnant:\nInfant",
                         TRUE ~ "Female, Pregnant:\nMaternal"),
         pop = factor(pop, levels = c("Male","Female, Non-pregnant", "Female, Pregnant:\nMaternal", "Female, Pregnant:\nInfant")),
         var_pop = case_when(var %in% c("f_still", "n_still") ~ "STILL",
                             var %in% c("f_lbw", "n_lbw") ~ "LBW",
                             var %in% c("f_pneu", "n_pneu") ~ "NP",
                             var %in% c("f_conj", "n_conj") ~ "ON",
                             var %in% c("f_cpp_tfi", "f_cpp_ep", "f_ep_tfi", "f_cpp_tfi_ep") ~ "Multiple",
                             TRUE ~ var)) |>
  mutate(var_pop = factor(var_pop, levels = c("f_ng","f_pid","f_cpp", "f_tfi", "f_ep", 
                                             "Multiple", "STILL", "LBW", "NP", "ON",
                                             "m_urethritis", "m_eds"),
                                   labels = c("NG","PID","CPP", "TFI", "EP", "Multiple",
                                             "STILL", "LBW", "NP", "ON", "UTS", "EDS"))) |>
  filter(!is.na(value)) |>
  group_by(lhs_num, sex, gest, pop, var_pop) |>
  summarise(value = sum(value)) |>
  mutate(prop = value/sum(value)) |>
  ungroup() 


p1 <- df_plot |>
  group_by(pop, var_pop) |>
  summarise(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) |>
  ggplot() +
  geom_col(aes(x = var_pop, y = median), width = 0.6) +
  facet_grid(~pop, scales = "free_x", space = "free") +
  mytheme +
  labs(x = "", y = "Total QALYs lost (in 1000)") +
  scale_y_continuous(labels = scales::label_number(scale = 1/1000)) +
  coord_cartesian(ylim = c(0, 65000))

p2 <- df_plot |>
  group_by(pop, var_pop) |>
  summarise(median = median(prop),
            lwr = quantile(prop, 0.025),
            upr = quantile(prop, 0.975)) |>
  ggplot() +
  geom_col(aes(x = var_pop, y = median), width = 0.6) +
  facet_grid(~pop, scales = "free_x", space = "free") +
  mytheme +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "", y = "Proportion QALYs lost") 

p1/p2 +
  plot_annotation(tag_level = "A") 

ggsave(".././plots/qaly_distribution.png", width = 19, height = 14, unit = "cm", dpi = 700)  

# qalys lost per case
qaly_seq_strata_2025 |> 
  group_by(lhs_num, test_num, sex, gest) |> 
  summarise(case = sum(case), annual_py = sum(annual_py)) |>
  left_join(df_plot) |>
  mutate(qaly_case = value/case*1000) |>
  group_by(pop, var_pop) |>
  summarise(median = median(qaly_case),
            lwr = quantile(qaly_case, 0.025),
            upr = quantile(qaly_case, 0.975)) |>
  ggplot() +
  geom_col(aes(x = var_pop, y = median)) +
  facet_grid(~pop, scales = "free_x", space = "free") +
  mytheme +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "", y = "QALYs lost per 1000 cases ") 

# qalys lost per py
qaly_seq_strata_2025 |> 
  group_by(lhs_num, test_num, sex, gest) |> 
  summarise(case = sum(case), annual_py = sum(annual_py)) |>
  left_join(df_plot) |>
  mutate(qaly_py = value/annual_py*1000) |>
  group_by(pop, var_pop) |>
  summarise(median = median(qaly_py),
            lwr = quantile(qaly_py, 0.025),
            upr = quantile(qaly_py, 0.975)) |>
  ggplot() +
  geom_col(aes(x = var_pop, y = median)) +
  facet_grid(~pop, scales = "free_x", space = "free") +
  mytheme +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "", y = "QALYs lost per 1000 person-years") 

# CALCS ----

# Percentage pregnancy sequalae due to still and lbw in infant

df_1 <- qaly_seq_strata_2025 |>
  filter(test_strat == "Screening", test_pop == "agyw_sa") |>
  pivot_longer(cols = c(starts_with("qaly")), names_to = "var") |>
  mutate(var = str_replace(var, "qaly_", "")) |>
  filter(!is.na(value)) |>
  group_by(lhs_num, sex, gest, var) |>
  summarise(value = sum(value)) 

df_1 |>
  filter(sex == "Female", gest == "Pregnant") |>
  left_join(df_1 |> 
              filter(sex == "Female", gest == "Pregnant") |>
              group_by(lhs_num) |> 
              summarise(tot = sum(value))) |>
  mutate(prop = value/tot) |>
  group_by(sex, gest, var) |>
  summarise(median = median(prop),
            lwr = quantile(prop, 0.025), 
            upr = quantile(prop, 0.975))


df_1 |>
  group_by(lhs_num) |>
  summarise(value = sum(value)) |>
  group_by() |>
  summarise(median = median(value),
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975))
