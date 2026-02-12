# INPUT ----

# UN WPP
pop <- read.csv("./parameters/pop_sag_calc.csv") |>
  mutate(age = factor(age_group, levels = c("15-24", "25-49"),
                      labels = c("15-24 y", "25-49 y")),
         sex = factor(sex, levels = c("Male","Female")))

# Model estimates
lhsnum <- read.csv("./posterior/post_lhs_nums.csv")$x

df_baseline <- data.table::fread("./posterior/post_baseline.csv") |>
  filter(lhs_num %in% lhsnum,
         test_strat == "h") |>
  mutate(sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
         age = factor(age, levels = c(1,2), labels = c("15-24 y", "25-49 y")),
         gest = factor(gest, levels = c(1,2), labels = c("Non-pregnant", "Pregnant")))

# PLOT ----

mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=0,r=0,b=0,l=0), "cm")),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        plot.tag.position = c(0.02, 1),
        axis.ticks = element_line(size = rel(1.0)))

pop |> 
  select(year, sex, age, pregnant, nonpregnant) |>
  pivot_longer(cols = c("pregnant", "nonpregnant"), names_to = "gest") |>
  mutate(gest = factor(gest, levels = c("nonpregnant", "pregnant"),
                       labels = c("Non-pregnant", "Pregnant"))) |>
  left_join(
    df_baseline |>
      group_by(lhs_num, year, sex, age, gest) |>
      summarise(N = sum(N), .groups = "drop") |>
      group_by(year, sex, age, gest) |>
      summarise(N = median(N))) |>
  pivot_longer(cols = c(value, N), names_to = "type") |>
  mutate(type = factor(type, levels = c("N", "value"), labels = c("Model Estimates","UN WPP Projections"))) |>
  ggplot(aes(x = year, y = value, colour= type, linetype=type)) +
  geom_line(size = 0.55) +
  facet_grid(gest~sex+age, scales = "free",
             labeller = label_wrap_gen(multi_line = FALSE)) +
  mytheme +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_colour_manual(values = c("#893957","grey40")) +
  scale_linetype_manual(values = c("solid", "solid")) +
  labs(x = "", y = "", colour = "", linetype = "")


ggsave(".././plots/demog.png", width = 18, height = 9, unit = "cm", dpi = 700)  


pop |> 
  select(year, sex, age, pop) |>
  group_by(year, sex) |>
  summarise(value = sum(pop)) |>
  left_join(
    df_baseline |>
      group_by(lhs_num, year, sex) |>
      summarise(N = sum(N), .groups = "drop") |>
      group_by(year, sex) |>
      summarise(N = median(N))) |>
  pivot_longer(cols = c(value, N), names_to = "type") |>
  mutate(type = factor(type, levels = c("value", "N"), labels = c("UN WPP", "Model"))) |>
  ggplot(aes(x = year, y = value, colour= type)) +
  geom_line(size = 1) +
  facet_grid(~sex, scales = "free", 
             labeller = label_wrap_gen(multi_line = FALSE)) +
  mytheme +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "", y = "", colour = "")
