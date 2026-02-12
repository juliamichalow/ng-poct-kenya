library(tidyverse)
library(patchwork)

# INPUTS ----

lhs_input <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "param")
prior_samples <- read.csv("./prior/prior_samples.csv", check.names = FALSE) |> mutate(lhs_num = row_number())
post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)
post_prev <- read.csv("./posterior/post_prev.csv", check.names = FALSE)
post_outcome <- read.csv("./posterior/post_outcome.csv", check.names = FALSE)

mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.2, "cm"),
        legend.position = "top",
        legend.margin = margin(unit(c(t=0,r=0,b=0,l=0), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.1), face="bold",
                                  margin = margin(c(t=1,r=0,b=1,l=0))),
        #strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))

order <- c("upsilon-3","upsilon-5","upsilon-4","upsilon-6",
           "q_low-1", "q_med_ratio-1", "q_high_ratio-1",
           
           "c_low-1","c_med_ratio-1","c_med_ratio-2","c_high_ratio-1","c_high_ratio-2","c_age_ratio-1","c_preg_ratio-1",
           
           "chi_low-1","chi_med_ratio-1","chi_high_ratio-1","chi_age_ratio-1",
           "chi_preg_ratio-1","chi_preg_ratio-2","chi_preg_ratio-3",
           
           "e-1",
           "gamma_low-1","gamma_sex_ratio-1","gamma_med_ratio-1","gamma_high_ratio-1","gamma_age_ratio-1","gamma_preg_ratio-1",
           
           "tau_m_low-1", "tau_m_sex_ratio-1",
           "n_l-1",
           "ap-1","ap-2",
           "zeta-1","zeta-2",
           
           "kappa-1","kappa-2",
           "phi-1","phi-2",
           "pi-1","sigma-1"
)

# PREP DATA ----

# columns to plot
plot_cols <- names(post_samples)[1:(which(names(post_samples) == "lhs_num") - 1)]

# x axis breaks
breaks <- post_samples |>
  select(all_of(c(plot_cols))) |>
  pivot_longer(cols = all_of(plot_cols), names_to = "parameter", values_to = "value") |>
  left_join(lhs_input |> select(var_strat, var_ggplot), 
            by = c("parameter"="var_strat")) |>
  group_by(var_ggplot) |>
  summarise(xmin = min(value),
            xmid = (max(value) + min(value))/2,
            xmax = max(value))

# full dataset
df_plot <- post_samples |>
  left_join(post_outcome, by = "lhs_num") |>
  pivot_longer(cols = all_of(plot_cols), names_to = "parameter", values_to = "value") |>
  left_join(lhs_input |> select(var_strat, var_ggplot), 
            by = c("parameter"="var_strat")) |>
  left_join(breaks, by = "var_ggplot") |>
  mutate(parameter = factor(parameter, levels = order),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])))

# PLOT ----

## Outcome stratified ----

# function
plot_distribution <- function(df_plot, outcome_name, title){
  
  colours1 <- c("#D32F2F", "#FFCCCB", "grey30", "#ACE5EE", "#33CCCC")
  #colours1 <- c("#33CCCC","#ACE5EE","grey30","#FFCCCB","#D32F2F")
  
  df_plot |>
    filter(outcome == outcome_name) |>
    ggplot() +
    #geom_density(aes(x = value, fill = outcome_quintile, colour = outcome_quintile), alpha = 0.3, linewidth = 0.3) +
    geom_density(aes(x = value, 
                     fill = forcats::fct_rev(outcome_quintile), 
                     colour = forcats::fct_rev(outcome_quintile)), 
                 alpha = 0.3, linewidth = 0.3) +
    facet_wrap(~var_ggplot, scales = "free", ncol = 7, labeller = label_parsed) +
    scale_fill_manual(name = title, values = colours1) +
    scale_colour_manual(name = title, values = colours1) +
    ggh4x::facetted_pos_scales(
      x = lapply(levels(df_plot$var_ggplot), function(v) {
        x <- df_plot[df_plot$var_ggplot == v, ][1,]  # Take first row for this var_ggplot
        scale_x_continuous(breaks = c(x$xmin, x$xmid, x$xmax),
                           labels = if(x$xmin < 10) {
                             sprintf("%.2f", c(x$xmin, x$xmid, x$xmax))
                           } else {
                             sprintf("%.1f", c(x$xmin, x$xmid, x$xmax))
                           })
      })) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    mytheme +
    theme(axis.text.x = element_text(angle = 90),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    # guides(fill = guide_legend(nrow = 1),
    #        colour = guide_legend(nrow = 1)) +
    guides(fill = guide_legend(nrow = 1, reverse = TRUE),
           colour = guide_legend(nrow = 1, reverse = TRUE)) +
    labs(x="",y="")
}

plot_distribution(df_plot, "prev_fm_overall", "Female prevalence")

plot_distribution(df_plot, "prev_fsw", "FSW prevalence")

plot_distribution(df_plot, "ratio_mf_prev", "Male-to-female prevalence ratio")

plot_distribution(df_plot, "ratio_mf_cases", "Male-to-female case ratio")

## Overall ----

df_post <- post_samples |> 
  select(all_of(c(plot_cols))) |>
  pivot_longer(cols = all_of(plot_cols), names_to = "parameter", values_to = "value") |>
  left_join(lhs_input |> select(var_strat, var_ggplot), 
            by = c("parameter"="var_strat")) |>
  left_join(breaks, by = "var_ggplot") |>
  mutate(parameter = factor(parameter, levels = order),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])))

df_prior <- prior_samples |>
  select(all_of(c(plot_cols))) |>
  pivot_longer(cols = all_of(plot_cols), names_to = "parameter", values_to = "value") |>
  left_join(lhs_input |> select(var_strat, var_ggplot), 
            by = c("parameter"="var_strat")) |>
  mutate(parameter = factor(parameter, levels = order),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])))

colours2 <- c("#616E61","#893957")

ggplot() +
  geom_density(data = df_prior, aes(x = value, fill = "Prior", colour = "Prior"), alpha = 0.3, linewidth = 0.3) +
  geom_density(data = df_post, aes(x = value, fill = "Posterior", colour = "Posterior"), alpha = 0.3, linewidth = 0.3) +
  facet_wrap(~var_ggplot, scales = "free", ncol = 7, labeller = label_parsed) +
  scale_fill_manual(name = "", breaks = c("Prior","Posterior"), values = colours2) +
  scale_colour_manual(name = "", breaks = c("Prior","Posterior"), values = colours2) +
  ggh4x::facetted_pos_scales(
    x = lapply(levels(df_post$var_ggplot), function(v) {
      x <- df_post[df_post$var_ggplot == v, ][1,]  # Take first row for this var_ggplot
      scale_x_continuous(breaks = c(x$xmin, x$xmid, x$xmax),
                         labels = if(x$xmin < 10) {
                           sprintf("%.2f", c(x$xmin, x$xmid, x$xmax))
                         } else {
                           sprintf("%.1f", c(x$xmin, x$xmid, x$xmax))
                         })
    })) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(fill = guide_legend(nrow = 1),
         colour = guide_legend(nrow = 1)) +
  labs(x="",y="")

ggsave(".././plots/posterior_all.png", width = 18, height = 17, unit = "cm", dpi = 400)
