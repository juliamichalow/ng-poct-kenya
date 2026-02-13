# NG prevalence in general population and fsw population

library(tidyverse)
library(patchwork)
library(glmmTMB)

df_fsw <- read.csv("./parameters/ngprev_fsw_adjusted.csv") |>
  select(study_id, study_group, year_mid, location_1, 
         pop, stratification, 
         age_group, age_range,
         hiv_status, hiv_prevalence,
         num, denom, specimencat, testcat, sens, spec,
         starts_with("adj_")) |>
  mutate(pop = "FSW",
         pop_cat = "FSW",
         sex = "Female",
         location = str_to_title(location_1),
         age_group = case_when(age_group == "adult" ~ "Adult", TRUE ~ age_group)) |>
  select(!c(location_1, adj_se))
    

df_gen <-  read.csv("./parameters/ngprev_gen_adjusted.csv") |>
  filter(sti == "NG") |>
  select(study_id, study_name, year_mid, location, 
         population, sex,
         age_group, age_range, hiv_status, num, denom, specimen, test, sens, spec,
         starts_with("adj_")) |>
  mutate(pop_cat = "Lower-risk",
         location = case_when(location == "Lumumba, Kisumu" ~ "Kisumu",
                              location == "Kisumu and surrounding villages" ~ "Kisumu",
                              location == "Siaya County" ~ "Siaya",
                              TRUE ~ location)) |>
  mutate(hiv_prevalence = case_when(hiv_status == "HIV negative" ~ 0,
                                    study_id == "Jespers 2014" & hiv_status == "Mixed" ~ NA,
                                    study_id == "Kerubo 2016" & hiv_status == "Mixed" ~ NA,
                                    study_id == "Masese 2017" & hiv_status == "Mixed" ~ NA,
                                    study_id == "Masha 2017" & hiv_status == "Mixed" ~ 13/202,
                                    study_id == "Mehta 2023" & hiv_status == "Mixed" ~ 7/436,
                                    study_id == "Oliver 2018" & hiv_status == "Mixed" ~ 67/457))

df_full <- full_join(df_fsw, df_gen,
                     by = c("study_id", "study_group" = "study_name","year_mid", "location", "pop" = "population", "pop_cat", "sex",
                            "age_group", "age_range","hiv_status", "hiv_prevalence",
                            "num", "denom",
                            "specimencat" = "specimen", "testcat" = "test", "sens", "spec",
                            "adj_num", "adj_denom", 
                            "adj_prev", "adj_prev_lwr", "adj_prev_upr")) |>
  mutate(hiv_status = factor(hiv_status, levels = c("HIV negative", "HIV positive","Mixed"),
                             labels = c("Negative","Positive","Non-stratified")),
         hiv_prevalence = as.numeric(hiv_prevalence),
         province = case_when(location %in% c("Western Kenya") ~ "Western",
                             location %in% c("Kisumu", "Siaya", "Siaya County", "Lumumba, Kisumu") ~ "Nyanza",
                             location %in% c("Mombasa", "Kilifi") ~ "Coast",
                             location %in% c("Thika") ~ "Central",
                             location %in% c("Nairobi") ~ "Nairobi",
                             location %in% c("Thika and Kisumu", "Nairobi and Mombasa") ~ "Multiple",
                             location %in% "NR" ~ "NR"),
         province = factor(province, levels = c("Coast", "Central", "Western", "Nyanza", "Nairobi", "Multiple", "NR")),
         # pop = case_when(pop %in% c("ANC attendees","FP attendees","GYN attendees") ~ "ANC/FP/GYN attendees", TRUE ~ pop),
         # pop = fct_relevel(pop, "ANC/FP/GYN attendees"),
         pop = factor(pop, 
                      levels = c("ANC attendees","FP attendees","GYN attendees","PHC/OPD attendees",
                                 "Students","Community members","HIV/STI prevention trial participants",
                                 "Population-representative survey participants","FSW"),
                      labels = c("ANC attendees","FP attendees","GYN attendees","PHC/OPD attendees",
                                 "Students","Community members","HIV/STI prevention trial participants",
                                 "Population-representative survey participants","Female sex workers")),
         group = fct_collapse(pop, 
                              "Overall" = c("Population-representative survey participants",
                                            "Community members", "PHC/OPD attendees"),
                              "Pregnant women" = c("ANC attendees","FP attendees","GYN attendees"),
                              "Higher-risk" = c("HIV/STI prevention trial participants")),
         group = factor(group, levels = c("Overall", "Pregnant women", "Higher-risk", "FSW")),
         risk = case_when(pop == "Female sex workers" ~ "FSW", TRUE ~ "Lower"),
         risk = fct_relevel(risk, "Lower"),
         sex = factor(sex, levels = c("Female", "Male")),
         age_group = factor(age_group, levels = c("Adult","Youth")),
         year_regression = year_mid - 2010)

# Leave out single male datapoint and rather use ratio to determine male prevalence
# df_full <- df_full |> filter(sex == "Female")

# Prepare table
df_analyse <- df_full |>
  filter(year_mid >= 2000) |>
  filter(!study_id == "Chanzu 2015") |> 
  mutate(prev = num/denom) 

df_analyse |>
  arrange(risk, sex, year_mid) |>
  select(sex, pop, year_mid, location, age_range, hiv_status, specimencat, testcat, sens, spec, 
         denom, prev, adj_prev, adj_prev_lwr, adj_prev_upr, study_id) |>
  # create latex row
  mutate(across(c(prev, adj_prev, adj_prev_lwr, adj_prev_upr), ~ sprintf("%.1f", . * 100)),
         across(c(sens, spec), ~ sprintf("%.1f",.)),
         adj_prev = paste0(adj_prev," (", adj_prev_lwr,"-",adj_prev_upr,")")) |>
  select(!c(adj_prev_lwr, adj_prev_upr)) |>
  rowwise() |>
  mutate(latex = paste0(paste(across(everything(), as.character), collapse = " & "), " \\\\")) # |> write.csv("../tables/table_studyprev.csv", row.names=FALSE)

# Prevalence per risk group ----

form2 <- cbind(adj_num,(adj_denom-adj_num)) ~ year_regression + risk + sex + age_group + hiv_status + (1 | study_id)

mod2 <- glmmTMB(form2, data = df_analyse, family = binomial(link="logit"))

sjPlot::tab_model(mod2)

df_pred <- expand.grid(sex = c("Female","Male"),
                       risk = c("FSW", "Lower"),
                       age_group = "Adult",
                       hiv_status = "Non-stratified",
                       year = c(2000:2025),
                       study_id = NA) |>
  mutate(year_regression = year - 2010) |>
  filter(!(sex == "Male" & risk %in% c("FSW")))

pred <- predict(mod2, newdata = df_pred, type="link", se.fit = TRUE, allow.new.levels=TRUE)

df_pred <- df_pred |>
  mutate(logit_prev = pred$fit, logit_se = pred$se.fit,
         prev = plogis(logit_prev),
         lwr = plogis(logit_prev - 1.96 * logit_se),
         upr = plogis(logit_prev + 1.96 * logit_se))

df_pred |>
  filter(year == 2020) |>
  mutate(est = paste0(round(prev*100,2), "% (", round(lwr*100,2), "-",round(upr*100,2),")")) |>
  select(risk, sex, hiv_status, age_group, year_regression, est)

# Plot
mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.2, "cm"),
        legend.position = "right",
        legend.margin = margin(unit(c(t=-5,r=5,b=5,l=5), "cm")),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold",
                                  margin = margin(unit(c(t=2,r=0,b=4,l=1), "cm")),
                                  hjust = 0, vjust = 0.5),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

df_pred <- df_pred |>
  mutate(plot_group = case_when(risk == "FSW" ~ "Female sex worker",
                                risk == "Lower" & sex == "Female" ~ "Female, overall",
                                risk == "Lower" & sex == "Male" ~ "Male, overall"),
         plot_group = factor(plot_group, levels = c("Female, overall", "Male, overall", "Female sex worker")))

df_analyse <- df_analyse |>
  mutate(plot_group = case_when(risk == "FSW" ~ "Female sex worker",
                                risk == "Lower" & sex == "Female" ~ "Female, overall",
                                risk == "Lower" & sex == "Male" ~ "Male, overall"),
         plot_group = factor(plot_group, levels = c("Female, overall", "Male, overall", "Female sex worker"))) 

df_analyse |>
  ggplot() +
  geom_ribbon(data = df_pred, aes(x = year, ymin = lwr, ymax= upr), alpha = 0.15) +
  geom_line(data = df_pred, aes(x = year, y = prev), size = 0.4) +
  geom_linerange(aes(x = year_mid, ymin = adj_prev_lwr, ymax = adj_prev_upr, colour = pop), 
                 size = 0.45) +
  geom_point(aes(x = year_mid, y = adj_prev, colour = pop, shape = hiv_status),
             size = 1, fill = "white") +
  facet_grid(~plot_group) +
  mytheme +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.margin = margin(unit(c(t=-5,r=5,b=5,l=5), "cm"))) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015, 2020, 2025)) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0,0.15)) +
  scale_shape_manual(values = c(21,19)) +
  scale_colour_manual(values = MetBrewer::met.brewer("Signac", n=9)) +
  guides(colour=guide_legend(order=2, nrow=3, title = "Population"),
         shape=guide_legend(order=1, title = "HIV status")) +
  labs(title = "", x = "", y = "")

ggsave("./parameters/ng_calibrate.png", width = 18, height = 9, unit = "cm", dpi = 700)  
