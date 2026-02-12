library(tidyverse)
library(patchwork)
library(ggtext)
library(ggh4x)
library(cowplot)


# INPUTS ----

impact_case <- read.csv("./intervention/impact_intervention_constrained_total.csv") |>
  filter(test_num == 25000) |>
  mutate(target_pop = factor(test_pop, levels = c("male_h","fsw","male_all","preg","agyw_sa"),
                             labels = c("CFSW","FSW","Men","Pregnant", "AGYW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))

impact_qaly <- read.csv("./qalys/impact_qaly_constrained_tot.csv") |>
  filter(test_num == 25000) |>
  mutate(target_pop = factor(test_pop, levels = c("male_h","fsw","male_all","preg","agyw_sa"),
                             labels = c("CFSW","FSW","Men","Pregnant", "AGYW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))

impact_case_max <- read.csv("./intervention/impact_intervention_unrestricted_total.csv") |>
  mutate(target_pop = factor(test_pop, levels = c("male_h","fsw","male_all","preg","agyw_sa"),
                             labels = c("CFSW","FSW","Men","Pregnant", "AGYW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))

impact_qaly_max <- read.csv("./qalys/impact_qaly_unrestricted_tot.csv") |>
  mutate(target_pop = factor(test_pop, levels = c("male_h","fsw","male_all","preg","agyw_sa"),
                             labels = c("CFSW","FSW","Men","Pregnant", "AGYW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))


mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=-10,r=10,b=0,l=10), "cm")),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(0.6)),
        panel.border = element_rect(color = "black", fill = NA, size = rel(0.9)))

# colours1 <- c("#9B2D5C","#D32F2F","#F15D52","#006D8F","#20A2AA")
# c("#B03060","#20B2AA")

col <- c("#DD6C11","#006D8F")

# PLOT 1 ----

plot_box <- function(data, yvar, ytitle, pos){
  
  # Calculate statistics for boxplot
  # Specific boxplot manually to use 95% CrI in whiskers rather than 1.5*IQR
  stats_df <- data |>
    group_by(target_pop, test_strat) |>
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
      aes(ymin = lwr,lower = q1, middle = median, upper = q3, ymax = upr, fill = test_strat),
      width = 0.6, fatten = 1, size = 0.25, alpha = 0.9,
      position = position_dodge(0.72)) +
    geom_label(
      #aes(label = paste0(sprintf("%.2f",median*100)), y = ypos, group = test_strat),
      #aes(label = ifelse(median < 0.002, paste0(sprintf("%.2f",median*100)), paste0(sprintf("%.1f",median*100))), y = ypos, group = test_strat),
      aes(label = paste0(sprintf("%.1f",median*100),"%"), y = ypos, group = test_strat),
      size = 7/.pt, label.size = NA, 
      label.padding = unit(0.02,"lines"),
      fill = "white", colour = "grey30",
      alpha = 0.75, hjust = 0.5, vjust = 0.5,
      position = position_dodge(0.72)) +
    scale_fill_manual(values = col) +
    mytheme +
    theme(legend.box = "vertical") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    #labs(x = "", y = ytitle, fill = "")
    labs(x = "", y = "", title = ytitle, fill = "")
}

p1 <- plot_box(impact_case |>
           mutate(y = cases_averted / cases_baseline),
         yvar = y,
         ytitle = "Percentage of infections averted", 
         pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "A")

p2 <- plot_box(impact_qaly |>
                 mutate(y = qaly_gained / qaly_baseline),
               yvar = y,
               ytitle = "Percentage of QALYs gained", 
               pos = 0.0021) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(tag = "B")

p3 <- plot_box(impact_case_max |>
                 mutate(y = cases_averted / cases_baseline),
               yvar = y,
               ytitle = "Percentage of infections averted", 
               pos = 0.019) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 0.6))  +
  labs(tag = "A")

p4 <- plot_box(impact_qaly_max |>
                 mutate(y = qaly_gained / qaly_baseline),
               yvar = y,
               ytitle = "Percentage of QALYs gained", 
               pos = 0.019) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 0.6))  +
  labs(tag = "B")


# PLOT 2 ----
# Number cases/QALYs losses averted per 100 test
plot_box <- function(data, yvar, ytitle, pos){
  
  # Calculate statistics for boxplot
  # Specific boxplot manually to use 95% CrI in whiskers rather than 1.5*IQR
  stats_df <- data |>
    group_by(target_pop, test_strat) |>
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
      aes(ymin = lwr,lower = q1, middle = median, upper = q3, ymax = upr, fill = test_strat),
      width = 0.6, fatten = 1, size = 0.25, alpha = 0.9,
      position = position_dodge(0.72)) +
    geom_label(
      aes(label = paste0(sprintf("%.1f",median)), y = ypos, group = test_strat),
      size = 7/.pt, label.size = NA, 
      label.padding = unit(0.02,"lines"),
      fill = "white", colour = "grey30",
      alpha = 0.75, hjust = 0.5, vjust = 0.5,
      position = position_dodge(0.72)) +
    scale_fill_manual(values = col) +
    mytheme +
    theme(legend.box = "vertical") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    #labs(x = "", y = ytitle, fill = "")
    labs(x = "", y = "", title = ytitle, fill = "")
}


p5 <- plot_box(impact_case |>
           mutate(y = cases_averted / test_used * 100),
         yvar = y,
         ytitle = "Infections averted per 100 tests", 
         pos = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 220)) +
  labs(tag = "C")

p6 <- plot_box(impact_qaly |>
           mutate(y = qaly_gained / test_used * 100),
         yvar = y,
         ytitle = "QALYs gained per 100 tests", 
         pos = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 22)) +
  labs(tag = "D")

p7 <- plot_box(impact_case_max |>
                 mutate(y = cases_averted / test_used * 100),
               yvar = y,
               ytitle = "Infections averted per 100 tests", 
               pos = 6.5) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 220)) +
  labs(tag = "C")

p8 <- plot_box(impact_qaly_max |>
                 mutate(y = qaly_gained / test_used * 100),
               yvar = y,
               ytitle = "QALYs gained per 100 tests", 
               pos = 0.65) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05))) +
  coord_cartesian(ylim = c(0, 22)) +
  labs(tag = "D")


# COMBINE

# Constrained POCT

top <- (p1 + p2) &
  theme(plot.title = element_text(size = 7.5*1.3, face = "bold"))

bot <- (p5 + p6) &
  theme(plot.title = element_text(size = 7.5*1.3, face = "bold"))


final_plot_1 <- (top / bot) +  
  plot_layout(guides = "collect") +
  plot_annotation(title = expression(underline(bold("Constrained POCT availability"))),
                  theme = theme(plot.title = element_text(size = 7.5*1.5, vjust = 2))) &
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

final_plot_1


ggsave(".././plots/intervention_constrained.png", width = 18, height = 17, unit = "cm", dpi = 700)  



# Unrestricted POCT

top <- (p3 + p4) &
  theme(plot.title = element_text(size = 7.5*1.3, face = "bold"))

bot <- (p7 + p8) &
  theme(plot.title = element_text(size = 7.5*1.3, face = "bold"))


final_plot_2 <- (top / bot) +  
  plot_layout(guides = "collect") +
  plot_annotation(title = expression(underline(bold("Unrestricted POCT availability"))),
                  theme = theme(plot.title = element_text(size = 7.5*1.5, vjust = 2))) &
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

final_plot_2

ggsave(".././plots/intervention_unrestricted.png", width = 18, height = 17, unit = "cm", dpi = 700)  


# PLOT 3 ----
# Number cases/QALYs losses averted per 100 test

case_eff <- 
  left_join(impact_case |>
              mutate(lim = cases_averted / test_used * 100) |>
              select(lhs_num, test_strat, target_pop, lim),
            
            impact_case_max |>
              mutate(max = cases_averted / test_used * 100) |>
              select(lhs_num, test_strat, target_pop, max)) |>
  pivot_longer(cols = c(lim, max), names_to = "test_group")

qaly_eff <- 
  left_join(impact_qaly |>
              mutate(lim = qaly_gained / test_used * 100) |>
              select(lhs_num, test_strat, target_pop, lim),
            
            impact_qaly_max |>
              mutate(max = qaly_gained / test_used * 100) |>
              select(lhs_num, test_strat, target_pop, max)) |>
  pivot_longer(cols = c(lim, max), names_to = "test_group")

plot_box <- function(data, yvar, ytitle, pos){
  
  # Calculate statistics for boxplot
  # Specific boxplot manually to use 95% CrI in whiskers rather than 1.5*IQR
  stats_df <- data |>
    group_by(target_pop, test_strat, test_group) |>
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
      aes(ymin = lwr,lower = q1, middle = median, upper = q3, ymax = upr, fill = test_group),
      width = 0.6, fatten = 1, size = 0.25, alpha = 0.9,
      position = position_dodge(0.72)) +
    geom_label(
      aes(label = paste0(sprintf("%.1f",median)), y = ypos, group = test_group),
      size = 7/.pt, label.size = NA, 
      label.padding = unit(0.02,"lines"),
      fill = "white", colour = "grey30",
      alpha = 0.75, hjust = 0.5, vjust = 0.5,
      position = position_dodge(0.72)) +
    facet_wrap(~test_strat) +
    scale_fill_manual(values = col, breaks = c("lim","max"), labels = c("Limited availability", "Maximum availability")) +
    mytheme +
    theme(legend.box = "vertical",
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(x = "", y = ytitle, fill = "")
}

plot_box(case_eff, value, "Infections averted per 100 test", pos = 5)
plot_box(qaly_eff, value, "QALYs gained per 100 test", pos = 0.5)


