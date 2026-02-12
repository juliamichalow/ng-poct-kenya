# PRCC between parameters and proportion of cases averted

library(qs)
library(sensitivity)
library(tidyverse)
library(patchwork)

# INPUTS ----

lhs_input <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "param")

post_impact <- read.csv("./intervention/impact_intervention_constrained_total.csv") |>
  filter(test_num == 25000)

post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)

post_prev <- read.csv("./posterior/post_prev.csv", check.names = FALSE)

post_outcome <- read.csv("./posterior/post_outcome.csv", check.names = FALSE) |>
  select(lhs_num, outcome, outcome_value) |>
  pivot_wider(names_from = outcome, values_from = outcome_value)


mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=0,r=0,b=5,l=0), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.text.y = element_text(size = rel(0.8)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        plot.tag.position = c(0.02, 1),
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

param_list <- names(post_samples)[1:(which(names(post_samples) == "lhs_num") - 1)]

df_param_outcome <- post_samples |>
  left_join(post_outcome, by = "lhs_num") 

df_param_impact <- post_samples |>
  left_join(post_impact, by = "lhs_num") 

# PRCC ----

calc_prcc <- function(df, params, outcome, n_bootstrap = 100) {
  
  # Create matrix X with parameter columns
  X <- as.data.frame(df[, params])
  
  # Rename columns to avoid problematic characters
  names(X) <- gsub("-", "_", params)
  
  # Create outcome vector
  y <- df[[outcome]]
  
  # Calculate PRCC
  prcc_results <- sensitivity::pcc(X = X, y = y, rank = TRUE, nboot = n_bootstrap)
  
  # Calculate PCC
  # pcc_results <- sensitivity::pcc(X = X, y = y, rank = FALSE, nboot = n_bootstrap)
  
  # Use PRCC when:
  # Unsure about linearity
  # Data might be skewed
  # You want to capture any monotonic relationship
  # Concerned about outliers
  
  # Use PCC when:
  # Confident relationships are linear
  # Data is normally distributed
  # You want to capture strictly linear effects
  
  prcc_results$PRCC |>
    rownames_to_column("parameter") |>
    rename(prcc = original,
           prcc_lwr = `min. c.i.`,
           prcc_upr = `max. c.i.`) |>
    mutate(sig = case_when((prcc_lwr >=0 & prcc_upr >=0) | 
                             (prcc_lwr <=0 & prcc_upr <=0) ~ "TRUE",
                           TRUE ~ "FALSE"),
           outcome= outcome)
}

## calibration targets ----

prcc_fm <- calc_prcc(df = df_param_outcome |> select(prev_fm_overall, all_of(param_list)),
                     params = param_list,
                     outcome = "prev_fm_overall")

prcc_fsw <- calc_prcc(df = df_param_outcome |> select(prev_fsw, all_of(param_list)),
                     params = param_list,
                     outcome = "prev_fsw")

prcc_ratio_prev <- calc_prcc(df = df_param_outcome |> select(ratio_mf_prev, all_of(param_list)),
                     params = param_list,
                     outcome = "ratio_mf_prev")

prcc_ratio_case <- calc_prcc(df = df_param_outcome |> select(ratio_mf_cases, all_of(param_list)),
                      params = param_list,
                      outcome = "ratio_mf_cases")

prcc_outcomes <- rbind(prcc_fm, prcc_fsw, prcc_ratio_prev, prcc_ratio_case) |>
  left_join(lhs_input |> select(var_strat, var_ggplot) |> mutate(var_strat = gsub("-", "_", var_strat)), 
            by = c("parameter"="var_strat")) |>
  mutate(parameter = factor(parameter, levels = gsub("-", "_", order)),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])),
         outcome = factor(outcome, 
                          levels = c("prev_fm_overall", "prev_fsw", "ratio_mf_prev", "ratio_mf_cases"),
                          labels = c("Female prevalence", "FSW prevalence", "Male-to-female prevalence ratio", "Male-to-female case ratio")))

write.csv(prcc_outcomes, "./sensitivity/prcc_outcomes.csv", row.names = FALSE)

## intervention impact ----
# loop through all the scenarios to calculate PRCC

# Define scenarios
test_populations <- c("preg", "agyw_sa", "fsw", "male_all", "male_h")
test_strategies <- c("Screening", "Diagnostic testing")
test_amounts <- c(25000)

# Initialize empty list to store results
all_results <- list()
i <- 1

# Loop through all scenarios
for (pop in test_populations) {
  for (strat in test_strategies) {
    for (amount in test_amounts) {
      # Filter data for current scenario
      df_scenario <- df_param_impact |>
        filter(test_pop == pop, 
               test_strat == strat, 
               test_num == amount)
      
      # Calculate PRCC for this scenario
      result_1 <- calc_prcc(df_scenario, params = param_list, "cases_averted")
      result_2 <- calc_prcc(df_scenario, params = param_list, "prop_cases_averted")
      
      # Add scenario information
      result_1 <- result_1 |>
        mutate(test_pop = pop,
               test_strat = strat,
               test_num = amount) |>
        relocate(outcome, test_strat, test_num, test_pop)
      
      result_2 <- result_2 |>
        mutate(test_pop = pop,
               test_strat = strat,
               test_num = amount) |>
        relocate(outcome, test_strat, test_num, test_pop)
      
      result <- rbind(result_1, result_2)
      
      # Store result in list
      all_results[[i]] <- result
      i <- i + 1
    }
  }
}

prcc_impact <- do.call(rbind, all_results) |>
  left_join(lhs_input |> select(var_strat, var_ggplot) |> mutate(var_strat = gsub("-", "_", var_strat)), 
            by = c("parameter"="var_strat")) |>
  mutate(parameter = factor(parameter, levels = gsub("-", "_", order)),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])),
         outcome = factor(outcome, 
                          levels = c("cases_averted", "prop_cases_averted"),
                          labels = c("Cases averted", "Proportion cases averted")),
         target_pop = factor(test_pop, levels = c("preg","agyw_sa","fsw","male_all","male_h"),
                             labels = c("PREG","AGYW","FSW","MEN", "CFSW")))

write.csv(prcc_impact, "./sensitivity/prcc_impact.csv", row.names = FALSE)

# PLOT ----

## Combined ----

colours <- c("grey70", "#B8E5CC", "#48B87C", "#2A8B65")

prcc_impact <- read.csv("./sensitivity/prcc_impact.csv") |>
  mutate(parameter = factor(parameter, levels = gsub("-", "_", order)),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])),
         target_pop = factor(test_pop, levels = c("preg","agyw_sa","fsw","male_all","male_h"),
                             labels = c("PREG","AGYW","FSW","MEN", "CFSW"))) |>
  mutate(cat = case_when(abs(prcc) >= 0.6 ~ "Strong",
                         abs(prcc) >= 0.4 ~ "Moderate",
                         abs(prcc) >= 0.2 ~ "Weak",
                         TRUE ~ "Negligible"),
         cat = factor(cat, levels = c("Negligible", "Weak", "Moderate", "Strong")))
  

# Diagnosis and screening combined
p1 <- prcc_impact |>
  filter(test_strat == "Diagnostic testing", outcome == "Proportion cases averted") |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW", "FSW", "MEN", "PREG", "AGYW"),
                             labels = c("CFSW", "FSW", "Men", "Pregnant", "AGYW"))) |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~target_pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6")) +
  labs(x = "", y = "", tag = "A", fill = "")

p2 <- prcc_impact |>
  filter(test_strat == "Screening", outcome == "Proportion cases averted") |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW", "FSW", "MEN", "PREG", "AGYW"),
                             labels = c("CFSW", "FSW", "Men", "Pregnant", "AGYW"))) |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~target_pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6")) +
  labs(x = "", y = "", tag = "B", fill = "")

p1 / p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(unit(c(t=-10,r=0,b=5,l=-5), "cm")),
        plot.margin = margin(unit(c(t=5,r=2,l=0,b=0),"cm")))

ggsave(".././plots/prcc_propcaseavert.png", width = 18, height = 25, unit = "cm", dpi = 700)  

## Separate ----

mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=0,r=0,b=5,l=0), "cm")),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = rel(1.2), face = "bold"),
        axis.text = element_text(size = rel(1.1)),
        axis.text.y = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1), face="bold"),
        legend.title = element_text(size = rel(1.1), face = "bold"),
        legend.text = element_text(size = rel(1.1)),
        strip.text = element_text(color="black", size = rel(1.3), face="bold"),
        strip.background = element_rect(color = NA, fill = NA),
        plot.tag = element_text(size=rel(1.4), face="bold"),
        axis.ticks = element_line(size = rel(1.0)))


p1 <- prcc_impact |>
  filter(test_strat == "Diagnostic testing", outcome == "Proportion cases averted") |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW", "FSW", "MEN", "PREG", "AGYW"),
                             labels = c("CFSW", "FSW", "Men", "Pregnant", "AGYW"))) |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~target_pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6")) +
  labs(x = "", y = "", tag = "A", fill = "")

p2 <- prcc_impact |>
  filter(test_strat == "Screening", outcome == "Proportion cases averted") |>
  mutate(target_pop = factor(target_pop, levels = c("CFSW", "FSW", "MEN", "PREG", "AGYW"),
                             labels = c("CFSW", "FSW", "Men", "Pregnant", "AGYW"))) |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~target_pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6")) +
  labs(x = "", y = "", tag = "B", fill = "")

ggsave(".././plots/prcc_propcaseavert_diag.png", p1, width = 18, height = 18, unit = "cm", dpi = 700)  
ggsave(".././plots/prcc_propcaseavert_screen.png", p2, width = 18, height = 18, unit = "cm", dpi = 700)  


# Analyse
prcc_impact |>
  filter(outcome == "Proportion cases averted") |>
  filter(cat %in% c("Moderate", "Strong")) |>
  arrange(test_strat, test_pop) |>
  select(test_strat, target_pop, parameter, var_ggplot, prcc, prcc_lwr, prcc_upr, sig, cat) 

