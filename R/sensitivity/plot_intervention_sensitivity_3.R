# Check how correlated parameters impact population ranking

post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)

impact_case <- read.csv("./intervention/impact_intervention_constrained_total.csv") |>
  filter(test_num == 25000) |>
  mutate(target_pop = factor(test_pop, levels = c("male_h","fsw","male_all","preg","agyw_sa"),
                             labels = c("CFSW","FSW","Men","Pregnant", "AGYW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening"))) |>
  left_join(post_samples)

impact_qaly <- read.csv("./qalys/impact_qaly_constrained_tot.csv") |>
  filter(test_num == 25000) |>
  mutate(target_pop = factor(test_pop, levels = c("male_h","fsw","male_all","preg","agyw_sa"),
                             labels = c("CFSW","FSW","Men","Pregnant", "AGYW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening"))) |>
  left_join(post_samples)

mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=-10,r=10,b=0,l=10), "cm")),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.title = element_text(size = rel(1.2), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

plot_diag <- function(data, yvar, ytitle, param, pos){
  
  # Calculate statistics for boxplot
  # Specific boxplot manually to use 95% CrI in whiskers rather than 1.5*IQR
  stats_df <- data |>
    filter(test_strat == "Diagnostic testing") |>
    group_by(target_pop, test_strat, param) |>
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
      aes(ymin = lwr,lower = q1, middle = median, upper = q3, ymax = upr, fill = param),
      width = 0.6, fatten = 1, size = 0.25, alpha = 0.9,
      position = position_dodge(0.72)) +
    geom_label(
      aes(label = paste0(sprintf("%.1f",median*100)), y = ypos, group = param),
      size = 7/.pt, label.size = NA, 
      label.padding = unit(0.02,"lines"),
      fill = "white", colour = "grey30",
      alpha = 0.75, hjust = 0.5, vjust = 0.5,
      position = position_dodge(0.72)) +
    mytheme +
    theme(legend.box = "vertical") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(x = "", y = ytitle, fill = "")
}

p1 <- plot_diag(impact_case |>
            mutate(y = cases_averted / cases_baseline,
                   param = cut(`zeta-1`, breaks = quantile(`zeta-1`, probs = seq(0, 1, 0.25)), 
                               include.lowest = TRUE,
                               labels = c("Q1", "Q2", "Q3", "Q4"))) |>
              filter(param %in% c("Q1", "Q4")),
          yvar = y,
          ytitle = "Percentage of infections averted",
          pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "A", fill = expression("Quartile of " * zeta[m] * " "))

p2 <- plot_diag(impact_qaly |>
            mutate(y = qaly_gained / qaly_baseline,
                   param = cut(`zeta-1`, breaks = quantile(`zeta-1`, probs = seq(0, 1, 0.25)), 
                               include.lowest = TRUE,
                               labels = c("Q1", "Q2", "Q3", "Q4"))) |>
              filter(param %in% c("Q1", "Q4")),
          yvar = y,
          ytitle = "Percentage of QALYs gained",
          pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "B", fill = expression("Quartile of " * zeta[m] * " "))

quantile(impact_case$`zeta-1`, probs = seq(0, 1, 0.25))

p3 <- plot_diag(impact_case |>
                  mutate(y = cases_averted / cases_baseline,
                         param = cut(`c_high_ratio-2`, breaks = quantile(`c_high_ratio-2`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of infections averted",
                pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "A", fill = expression("Quartile of " * rc[fh:fi] * " "))

p4 <- plot_diag(impact_qaly |>
                  mutate(y = qaly_gained / qaly_baseline,
                         param = cut(`c_high_ratio-2`, breaks = quantile(`c_high_ratio-2`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of QALYs gained",
                pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "B", fill = expression("Quartile of " * rc[fh:fi] * " "))

p5 <- plot_diag(impact_case |>
                  mutate(y = cases_averted / cases_baseline,
                         param = cut(`ap-2`, breaks = quantile(`ap-2`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of infections averted",
                pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "A", fill = expression("Quartile of " * p[f] * " "))

p6 <- plot_diag(impact_qaly |>
                  mutate(y = qaly_gained / qaly_baseline,
                         param = cut(`ap-2`, breaks = quantile(`ap-2`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of QALYs gained",
                pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "B", fill = expression("Quartile of " * p[f] * " "))

p7 <- plot_diag(impact_case |>
                  mutate(y = cases_averted / cases_baseline,
                         param = cut(`kappa-2`, breaks = quantile(`kappa-2`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of infections averted",
                pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "A", fill = expression("Quartile of " * kappa[f] * " "))

p8 <- plot_diag(impact_qaly |>
                  mutate(y = qaly_gained / qaly_baseline,
                         param = cut(`kappa-2`, breaks = quantile(`kappa-2`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of QALYs gained",
                pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "B", fill = expression("Quartile of " * kappa[f] * " "))

p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p3 + p4 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p5 + p6 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p7 + p8 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


plot_screen <- function(data, yvar, ytitle, param, pos){
  
  # Calculate statistics for boxplot
  # Specific boxplot manually to use 95% CrI in whiskers rather than 1.5*IQR
  stats_df <- data |>
    filter(test_strat == "Screening") |>
    group_by(target_pop, test_strat, param) |>
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
      aes(ymin = lwr,lower = q1, middle = median, upper = q3, ymax = upr, fill = param),
      width = 0.6, fatten = 1, size = 0.25, alpha = 0.9,
      position = position_dodge(0.72)) +
    geom_label(
      aes(label = paste0(sprintf("%.1f",median*100)), y = ypos, group = param),
      size = 7/.pt, label.size = NA, 
      label.padding = unit(0.02,"lines"),
      fill = "white", colour = "grey30",
      alpha = 0.75, hjust = 0.5, vjust = 0.5,
      position = position_dodge(0.72)) +
    mytheme +
    theme(legend.box = "vertical") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(x = "", y = ytitle, fill = "")
}

p1 <- plot_screen(impact_case |>
                  mutate(y = cases_averted / cases_baseline,
                         param = cut(`upsilon-5`, breaks = quantile(`upsilon-5`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of infections averted",
                pos = 0.001) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.02)) +
  labs(tag = "A", fill = expression("Quartile of " * upsilon[mh] * " "))

p2 <- plot_screen(impact_qaly |>
                  mutate(y = qaly_gained / qaly_baseline,
                         param = cut(`upsilon-5`, breaks = quantile(`upsilon-5`, probs = seq(0, 1, 0.25)), 
                                     include.lowest = TRUE,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                  filter(param %in% c("Q1", "Q4")),
                yvar = y,
                ytitle = "Percentage of QALYs gained",
                pos = 0.001) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.02)) +
  labs(tag = "B", fill = expression("Quartile of " * upsilon[mh] * " "))


p3 <- plot_screen(impact_case |>
                    mutate(y = cases_averted / cases_baseline,
                           param = cut(`upsilon-6`, breaks = quantile(`upsilon-6`, probs = seq(0, 1, 0.25)), 
                                       include.lowest = TRUE,
                                       labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                    filter(param %in% c("Q1", "Q4")),
                  yvar = y,
                  ytitle = "Percentage of infections averted",
                  pos = 0.001) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.02)) +
  labs(tag = "A", fill = expression("Quartile of " * upsilon[fh] * " "))

p4 <- plot_screen(impact_qaly |>
                    mutate(y = qaly_gained / qaly_baseline,
                           param = cut(`upsilon-6`, breaks = quantile(`upsilon-6`, probs = seq(0, 1, 0.25)), 
                                       include.lowest = TRUE,
                                       labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                    filter(param %in% c("Q1", "Q4")),
                  yvar = y,
                  ytitle = "Percentage of QALYs gained",
                  pos = 0.001) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.02)) +
  labs(tag = "B", fill = expression("Quartile of " * upsilon[fh] * " "))



p1 <- plot_screen(impact_case |>
                    mutate(y = cases_averted / cases_baseline,
                           param = cut(`c_high_ratio-2`, breaks = quantile(`c_high_ratio-2`, probs = seq(0, 1, 0.25)), 
                                       include.lowest = TRUE,
                                       labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                    filter(param %in% c("Q1", "Q4")),
                  yvar = y,
                  ytitle = "Percentage of infections averted",
                  pos = 0.001) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.02)) +
  labs(tag = "A", fill = expression("Quartile of " * rc[fh:fi] * " "))

p2 <- plot_screen(impact_qaly |>
                    mutate(y = qaly_gained / qaly_baseline,
                           param = cut(`c_high_ratio-2`, breaks = quantile(`c_high_ratio-2`, probs = seq(0, 1, 0.25)), 
                                       include.lowest = TRUE,
                                       labels = c("Q1", "Q2", "Q3", "Q4"))) |>
                    filter(param %in% c("Q1", "Q4")),
                  yvar = y,
                  ytitle = "Percentage of QALYs gained",
                  pos = 0.001) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(0, 0.02)) +
  labs(tag = "B", fill = expression("Quartile of " * rc[fh:fi] * " "))

p1 + p2