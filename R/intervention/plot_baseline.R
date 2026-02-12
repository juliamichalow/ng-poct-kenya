## PREP DATA ----

lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

intervention_strata <- data.table::fread("./posterior/post_intervention_constrained_strata.csv") |>
  filter(lhs_num %in% lhsnum)

# Percentage of tests among adult population
intervention_strata |> 
  filter(test_num == 0, test_pop == "fsw", test_strat == "h", year >= 2025) |>
  group_by(lhs_num, year) |> 
  summarise(N = sum(N), NS = sum(NS)) |>
  mutate(prop_N = 25000/N,
         prop_NS = 25000/NS) |>
  ungroup() |>
  group_by(year) |>
  summarise(prop_N = median(prop_N)*100,
            prop_NS = median(prop_NS)*100)



inf_strata <- intervention_strata |>
  filter(test_num == 0, test_pop == "fsw", test_strat == "h",
         year == 2025) |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
         risk = factor(risk, levels = c(1,2,3), labels = c("Low", "Medium", "High")),
         age = factor(age, levels = c(1,2), labels = c("Young", "Old")),
         gest = factor(gest, levels = c(1,2), labels = c("Non-pregnant", "Pregnant"))) |>
  filter(!(sex == "Male" & gest == "Pregnant"))

rm(intervention_strata)

qaly_strata <- read.csv("./qalys/impact_qaly_constrained_tot_strata_2025.csv") |>
  filter(test_pop == "agyw_sa",
         test_strat == "Screening") |>
  mutate(sex = factor(sex), risk = factor(risk), age = factor(age), gest = factor(gest)) |>
  select(!c(test_pop, test_strat, test_num, test_used)) |>
  left_join(inf_strata |> select(lhs_num,sex,risk,age,gest,NS))



# Prep data per population group
# Prevalence
calc_prev <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(I)/sum(NS))
}

df_prev <- rbind(
  calc_prev(inf_strata, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_prev(inf_strata, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "Pregnant"),
  calc_prev(inf_strata, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_prev(inf_strata, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "Men"),
  calc_prev(inf_strata, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW","FSW","Men","Pregnant","AGYW")))

df_prev |>
  group_by(target_pop) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 

# incident cases
calc_case <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(n_I))
}

df_case <- rbind(
  calc_case(inf_strata, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_case(inf_strata, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "Pregnant"),
  calc_case(inf_strata, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_case(inf_strata, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "Men"),
  calc_case(inf_strata, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW","FSW","Men","Pregnant","AGYW")))

# Incidence per 100 population
calc_inc <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(n_I)/sum(NS)*100)
}

df_inc <- rbind(
  calc_inc(inf_strata, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_inc(inf_strata, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "Pregnant"),
  calc_inc(inf_strata, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_inc(inf_strata, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "Men"),
  calc_inc(inf_strata, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW","FSW","Men","Pregnant","AGYW")))

df_inc |>
  group_by(target_pop) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 


# Incident cases as % of total 
calc_propinc <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(case = sum(n_I)) |>
    left_join(data |>
                group_by(lhs_num) |>
                summarise(tot = sum(n_I))) |>
    mutate(value = case/tot)
}

df_propinc <- rbind(
  calc_propinc(inf_strata, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_propinc(inf_strata, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "Pregnant"),
  calc_propinc(inf_strata, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_propinc(inf_strata, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "Men"),
  calc_propinc(inf_strata, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW","FSW","Men","Pregnant","AGYW")))

df_propinc |>
  group_by(target_pop) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 


# QALYs lost per 1000 population
calc_qaly <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(value = sum(qaly_baseline)/sum(NS)*1000)
}

df_qaly <- rbind(
  calc_qaly(qaly_strata, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_qaly(qaly_strata, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "Pregnant"),
  calc_qaly(qaly_strata, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_qaly(qaly_strata, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "Men"),
  calc_qaly(qaly_strata, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW","FSW","Men","Pregnant","AGYW")))

df_qaly |>
  group_by(target_pop) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 

# QALYs lost as % of total 
calc_propqaly <- function(data, ...) {
  
  data |>
    filter(...) |>
    group_by(lhs_num) |>
    summarise(qaly = sum(qaly_baseline)) |>
    left_join(data |>
                group_by(lhs_num) |>
                summarise(tot = sum(qaly_baseline))) |>
    mutate(value = qaly/tot)
}

df_propqaly <- rbind(
  calc_propqaly(qaly_strata, sex == "Female", age == "Young", gest == "Non-pregnant") |> mutate(target_pop = "AGYW"),
  calc_propqaly(qaly_strata, sex == "Female", gest == "Pregnant") |> mutate(target_pop = "Pregnant"),
  calc_propqaly(qaly_strata, sex == "Female", risk == "High") |> mutate(target_pop = "FSW"),
  calc_propqaly(qaly_strata, sex == "Male", gest == "Non-pregnant") |> mutate(target_pop = "Men"),
  calc_propqaly(qaly_strata, sex == "Male", gest == "Non-pregnant", risk == "High") |> mutate(target_pop = "CFSW")) |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW","FSW","Men","Pregnant","AGYW")))

df_propqaly |>
  group_by(target_pop) |>
  summarise(median = median(value), 
            lwr = quantile(value, 0.025), 
            upr = quantile(value, 0.975)) 

# PLOT ----

mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=0,r=5,b=0,l=5), "cm")),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title.x = element_text(size = rel(1.1), face="bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(0.6)),
        panel.border = element_rect(color = "black", fill = NA, size = rel(0.9)))

# PLOT ----

label_prev <- df_prev |>
  group_by(target_pop) |>
  summarise(median = median(value),
            q3 = quantile(value, 0.75),
            ypos = q3 + 0.005)

p1 <- df_prev |>
  ggplot(aes(x = target_pop, y = value)) +
  geom_boxplot(width = 0.5, fatten = 1, size = 0.25, alpha = 0.7, outliers = FALSE, fill = "grey40") +
  mytheme +
  geom_label(data = label_prev, 
    aes(label = paste0(sprintf("%.1f",median*100)), y = ypos),
    size = 7.5/.pt, label.size = NA, 
    label.padding = unit(0.035,"lines"),
    fill = "white", colour = "grey30",
    alpha = 0.75, hjust = 0.5, vjust = 0.5) +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0, 0.1),
                     #breaks = c(0.001, 0.003, 0.01, 0.03, 0.1),
                     #trans = "log",
                     #breaks = scales::log_breaks(n = 5),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "", y = "Prevalence (%)", fill = "")

label_inc <- df_inc |>
  group_by(target_pop) |>
  summarise(median = median(value),
            q3 = quantile(value, 0.75)) |>
  mutate(ypos = ifelse(target_pop == "CFSW", q3 + 0.6, q3 + 1))

p2 <- df_inc |>
  ggplot(aes(x = target_pop, y = value)) +
  geom_boxplot(width = 0.5, fatten = 1, size = 0.25, alpha = 0.7, outliers = FALSE, fill = "grey40") +
  mytheme +
  geom_label(data = label_inc, 
             aes(label = paste0(sprintf("%.1f", median)), y = ypos),
             size = 7.5/.pt, label.size = NA, 
             label.padding = unit(0.035,"lines"),
             fill = "white", colour = "grey30",
             alpha = 0.75, hjust = 0.5, vjust = 0.5) +
  scale_y_continuous(limits = c(0, 25),
                     breaks = c(0, 5, 10, 15, 20, 25),
                     expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "", y = "", title = "Incident infections per 100 population", fill = "")


p2

label_qaly <- df_qaly |>
  group_by(target_pop) |>
  summarise(median = median(value),
            q3 = quantile(value, 0.75)) |>
  mutate(ypos = ifelse(target_pop == "Pregnant", q3 + 3, q3 + 5))

p3 <- df_qaly |>
  ggplot(aes(x = target_pop, y = value)) +
  geom_boxplot(width = 0.5, fatten = 1, size = 0.25, alpha = 0.7, outliers = FALSE, fill = "grey40") +
  mytheme +
  geom_label(data = label_qaly, 
             aes(label = paste0(sprintf("%.1f", median)), y = ypos),
             size = 7.5/.pt, label.size = NA, 
             label.padding = unit(0.035,"lines"),
             fill = "white", colour = "grey30",
             alpha = 0.75, hjust = 0.5, vjust = 0.5) +
  scale_y_continuous(limits = c(0, 150),
                     breaks = c(0, 30, 60, 90, 120, 150),
                     expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "", y = "", title = "QALYs lost per 1000 population", fill = "")


label_propinc <- df_propinc |>
  group_by(target_pop) |>
  summarise(median = median(value),
            q3 = quantile(value, 0.75),
            ypos = q3 + 0.04)

p4 <- df_propinc |>
  ggplot(aes(x = target_pop, y = value)) +
  geom_boxplot(width = 0.5, fatten = 1, size = 0.25, alpha = 0.7, outliers = FALSE, fill = "grey40") +
  mytheme +
  geom_label(data = label_propinc, 
             aes(label = paste0(sprintf("%.1f",median*100), "%"), y = ypos),
             size = 7.5/.pt, label.size = NA, 
             label.padding = unit(0.035,"lines"),
             fill = "white", colour = "grey30",
             alpha = 0.75, hjust = 0.5, vjust = 0.5) +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "", y = "", title = "Proportion of total incident infections", fill = "")


label_propqaly<- df_propqaly |>
  group_by(target_pop) |>
  summarise(median = median(value),
            q3 = quantile(value, 0.75),
            ypos = q3 + 0.04)

p5 <- df_propqaly |>
  ggplot(aes(x = target_pop, y = value)) +
  geom_boxplot(width = 0.5, fatten = 1, size = 0.25, alpha = 0.7, outliers = FALSE, fill = "grey40") +
  mytheme +
  geom_label(data = label_propqaly, 
             aes(label = paste0(sprintf("%.1f",median*100), "%"), y = ypos),
             size = 7.5/.pt, label.size = NA, 
             label.padding = unit(0.035,"lines"),
             fill = "white", colour = "grey30",
             alpha = 0.75, hjust = 0.5, vjust = 0.5) +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "", y = "", title = "Proportion of total QALYs lost", fill = "")



(p4 + p2) / (p5 + p3) +
  plot_annotation(tag_level = "A")


ggsave(".././plots/baseline.png", width = 18, height = 15.8, unit = "cm", dpi = 700)  

