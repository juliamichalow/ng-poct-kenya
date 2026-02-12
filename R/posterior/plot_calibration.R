library(tidyverse)
library(patchwork)

# INPUTS ----

# posterior outcomes
post_outcome <- read.csv("./posterior/post_outcome.csv", check.names = FALSE) |>
  filter(!outcome == "prev_m_overall") |>
  mutate(category = case_when(outcome %in% c("prev_fm_overall", "prev_fsw") ~ "Prevalence",
                              outcome %in% c("ratio_mf_prev", "ratio_mf_cases") ~ "Male-to-female ratio"),
         category = fct_relevel(category, "Prevalence"),
         outcome = factor(outcome, levels = c("prev_fm_overall", "prev_fsw", "ratio_mf_prev", "ratio_mf_cases"),
                          labels = c("Sexually active females", "Female sex workers", "Prevalence", "Symptomatic case rate")))

# calibration targets
targets <- data.frame(category = c("Prevalence", "Prevalence", "Male-to-female ratio", "Male-to-female ratio"),
                      outcome = c("Sexually active females", "Female sex workers", "Prevalence", "Symptomatic case rate"),
                      median = c(0.0225, 0.0425, 0.8, 0.475),
                      lwr = c(0.01, 0.02, 0.6, 0.35),
                      upr = c(0.035, 0.065, 1, 0.6)) |>
  mutate(category = fct_relevel(category, "Prevalence"))


mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.4, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=0,r=-5,b=-5,l=0), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold", vjust = 1),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

# PLOTS ----

# main version

# plot posterior outcomes against targets
p1 <- post_outcome |>
  group_by(category, outcome) |>
  summarise(median = median(outcome_value),
                lwr = quantile(outcome_value, 0.025),
                upr = quantile(outcome_value, 0.975)) |>
  mutate(type = "Posterior median (UI)") |>
  filter(category == "Prevalence") |>
  ggplot() +
  geom_violin(data = post_outcome |> filter(category == "Prevalence"),
              aes(x = outcome, y = outcome_value, fill = "Posterior distribution"), trim = TRUE,
              size = 0.35, alpha = 0.35, colour = NA) +
  geom_linerange(aes(x = outcome, ymin = lwr, ymax = upr, colour = type),
                 size = 0.45) +
  geom_point(aes(x = outcome, y = median, colour = type), shape = "diamond",
             size = 2) +
  ggh4x::facet_wrap2(~category, scales = "free", nrow = 1) +
  scale_fill_manual(values = "#893957") +
  scale_colour_manual(values = c("#893957")) + 
  mytheme +
  scale_y_continuous(limits = c(0, 0.08), labels = scales::label_percent()) +
  labs(x = "", y = "", fill = "", colour = "")

p2 <- post_outcome |>
  group_by(category, outcome) |>
  summarise(median = median(outcome_value),
            lwr = quantile(outcome_value, 0.025),
            upr = quantile(outcome_value, 0.975)) |>
  mutate(type = "Posterior median (UI)") |>
  filter(category == "Male-to-female ratio") |>
  ggplot() +
  geom_violin(data = post_outcome |> filter(category == "Male-to-female ratio"),
              aes(x = outcome, y = outcome_value, fill = "Posterior distribution"), trim = TRUE,
              size = 0.35, alpha = 0.35, colour = NA) +
  geom_linerange(aes(x = outcome, ymin = lwr, ymax = upr, colour = type),
                 size = 0.45) +
  geom_point(aes(x = outcome, y = median, colour = type), shape = "diamond",
             size = 2) +
  ggh4x::facet_wrap2(~category, scales = "free", nrow = 1) +
  scale_fill_manual(values = "#893957") +
  scale_colour_manual(values = c("#893957")) + 
  mytheme +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_comma()) +
  labs(x = "", y = "", fill = "", colour = "")

p1 + p2 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_level = "A") &
  theme(legend.position = "bottom")

ggsave(".././plots/posterior_target.png", width = 18, height = 9, unit = "cm", dpi = 700)  


# alternate?
targets |> mutate(type = "Calibration target") |>
  bind_rows(
    post_outcome |>
      group_by(category, outcome) |>
      summarise(median = median(outcome_value),
                lwr = quantile(outcome_value, 0.025),
                upr = quantile(outcome_value, 0.975)) |>
      mutate(type = "Posterior estimate")) |>
  ggplot() +
  geom_violin(data = post_outcome,
              aes(x = outcome, y = outcome_value, fill = "Posterior estimate distribution"), trim = TRUE,
              size = 0.35, alpha = 0.35, colour = "#893957") +
  geom_linerange(aes(x = outcome, ymin = lwr, ymax = upr, colour = type),
                 size = 0.5,
                 position = position_dodge(0.2)) +
  geom_point(aes(x = outcome, y = median, colour = type),
             size = 1.4,
             position = position_dodge(0.2)) +
  ggh4x::facet_wrap2(~category, scales = "free", nrow = 1) +
  scale_fill_manual(values = "#893957") +
  scale_colour_manual(values = c("grey20","#893957")) + 
  mytheme +
  ggh4x::facetted_pos_scales(
    y = list(
      category == "Prevalence" ~ scale_y_continuous(limits = c(0, 0.08), labels = scales::label_percent()),
      category == "Male-to-female ratio" ~ scale_y_continuous(limits = c(0, 1), labels = scales::label_comma()))) +
  labs(x = "", y = "", fill = "", colour = "")


