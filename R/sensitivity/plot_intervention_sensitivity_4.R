# Assess sensitivity of prioritisation results to:
# LBW duration and disutility of stillbirth

library(tidyverse)
library(patchwork)
library(ggtext)
library(ggh4x)

# INPUTS ----

impact_qaly_sens <- read.csv("./qalys/impact_qaly_constrained_tot_sensitivity.csv") |>
  mutate(target_pop = factor(test_pop, levels = c("male_h","fsw","male_all","preg","agyw_sa"),
                             labels = c("CFSW","FSW","Men","Pregnant", "AGYW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")),
         sens = factor(sens, levels = c("Primary", "Conservative", "No SB"))) |>
  mutate(group = paste(test_strat, sens),
         group = factor(group, levels = c(
           "Diagnostic testing Primary", 
           "Diagnostic testing Conservative", 
           "Diagnostic testing No SB",
           "Screening Primary", 
           "Screening Conservative", 
           "Screening No SB"
         )))


mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=0,r=10,b=0,l=10), "cm")),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

col <- c("#DD6C11","#DD6C11","#DD6C11","#006D8F","#006D8F","#006D8F")

# PLOT ----
## QALYs gained ----

plot_box <- function(data, yvar, ytitle, pos){
  
  # Calculate statistics for boxplot
  # Specific boxplot manually to use 95% CrI in whiskers rather than 1.5*IQR
  stats_df <- data |>
    group_by(target_pop, group, test_strat, sens) |>
    summarise(median = median({{yvar}}),
              q1 = quantile({{yvar}}, 0.25),
              q3 = quantile({{yvar}}, 0.75),
              lwr = quantile({{yvar}}, 0.025),  
              upr = quantile({{yvar}}, 0.975)) |>
    mutate(ypos = q3 + pos)
  
  stats_df |>
    ggplot(aes(x = target_pop)) +
    geom_boxplot(
      stat = "identity",
      aes(ymin = lwr,lower = q1, middle = median, upper = q3, ymax = upr, 
          fill = group, alpha = sens),
      width = 0.6, fatten = 1, size = 0.25,
      position = position_dodge(0.72)) +
    geom_label(
      aes(label = paste0(sprintf("%.1f",median*100)), y = ypos, group = group),
      size = 7/.pt, label.size = NA,
      label.padding = unit(0.02,"lines"),
      fill = "white", colour = "grey30",
      alpha = 0.75, hjust = 0.5, vjust = 0.5,
      position = position_dodge(0.72)) +
    mytheme +
    theme(legend.box = "vertical") +
    scale_fill_manual(values = col, guide = "none") +
    scale_alpha_manual(values = c(1, 0.5, 0.25), guide = "none") +  
    # new scale fill
    ggnewscale::new_scale_fill() +
    # Add invisible geoms for legend only
    # Two separate geom_point calls with different aesthetics
    geom_point(aes(y = -Inf, color = sens), shape = 22, size = 0) +
    geom_point(aes(y = -Inf, fill = test_strat), shape = 22, size = 0, stroke = 0, color = NA) +
    scale_colour_manual(values = c("grey10", "grey50", "grey70"),
                        breaks = c("Primary", "Conservative", "No SB"),
                        labels = c("Primary", "Conservative", "No stillbirth"),
                        guide = guide_legend(override.aes = list(shape = 15, size = 2.4),
                                             nrow = 1, order = 2, byrow = FALSE)) +
    scale_fill_manual(values = c("#DD6C11","#006D8F"),
                      breaks = c("Diagnostic testing", "Screening"),
                      guide = guide_legend(override.aes = list(shape = 22, size = 3),
                                           nrow = 1, order = 1, byrow = FALSE)) +
    labs(x = "", y = ytitle, fill = "Test strategy ", colour = "Parameters ") 
}



# All variations
plot_box(impact_qaly_sens |>
           mutate(y = qaly_gained / qaly_baseline),
         yvar = y,
         pos = 0.0016,
         ytitle = "Percentage of QALYs gained") +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.063))


ggsave(".././plots/qaly_impact_sensitivity.png", width = 18, height = 9, unit = "cm", dpi = 700)  


## Distribution of QALYs -----

df_plot <- data.table::fread("./qalys/impact_qaly_constrained_seq_strata_2025_sensitivity.csv") |>
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
  group_by(sens, lhs_num, sex, gest, pop, var_pop) |>
  summarise(value = sum(value)) |>
  mutate(prop = value/sum(value)) |>
  ungroup() 


p1 <- df_plot |>
  group_by(sens, pop, var_pop) |>
  summarise(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) |>
  ggplot() +
  geom_col(aes(x = var_pop, y = median, fill=sens, 
               group = sens),
               position = position_dodge(0.9), width = 0.6) +
  facet_grid(~pop, scales = "free_x", space = "free") +
  mytheme +
  labs(x = "", y = "Total QALYs lost (in 1000)") +
  scale_y_continuous(labels = scales::label_number(scale = 1/1000)) +
  coord_cartesian(ylim = c(0, 65000))

p2 <- df_plot |>
  group_by(sens, pop, var_pop) |>
  summarise(median = median(prop),
            lwr = quantile(prop, 0.025),
            upr = quantile(prop, 0.975)) |>
  ggplot() +
  geom_col(aes(x = var_pop, y = median, fill=sens), 
           position = position_dodge(0.9), width = 0.6) +
  facet_grid(~pop, scales = "free_x", space = "free") +
  mytheme +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "", y = "Proportion QALYs lost") 

p1/p2 +
  plot_annotation(tag_level = "A") 
