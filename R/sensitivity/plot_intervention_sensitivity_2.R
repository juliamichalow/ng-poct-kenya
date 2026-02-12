# PRCC between parameters and population ranking

# INPUTS ----

lhs_input <- readxl::read_excel("./parameters/parameters.xlsx", sheet = "param")

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


impact_case <- read.csv("./intervention/impact_intervention_constrained_total.csv") |>
  filter(test_num == 25000) |>
  mutate(target_pop = factor(test_pop, levels = c("preg","agyw_sa","fsw","male_all","male_h"),
                             labels = c("PREG","AGYW","FSW","MEN", "CFSW")),
         test_strat = factor(test_strat, levels = c("Diagnostic testing","Screening")))


post_samples <- read.csv("./posterior/post_samples.csv", check.names = FALSE)

mytheme <- theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(unit(c(t=0,r=0,b=0,l=0), "cm")),
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
        axis.ticks = element_line(size = rel(1.0)))

# CALC RANKS ----

# Add rankings for each parameter set and strategy
rank_populations <- function(AGYW, FSW, MEN, CFSW, PREG) {
  values <- c(AGYW, FSW, MEN, CFSW, PREG)
  rank(-values, ties.method = "min")
}

# Reshape impact data to wide format for easier ranking
impact_wide <- impact_case |> 
  select(lhs_num, test_strat, target_pop, prop_cases_averted) |>
  pivot_wider(id_cols = c(lhs_num, test_strat),
              names_from = target_pop,
              values_from = prop_cases_averted)

# Apply ranking
impact_ranked <- impact_wide %>%
  rowwise() %>%
  mutate(
    AGYW_rank = rank_populations(AGYW, FSW, MEN, CFSW, PREG)[1],
    FSW_rank = rank_populations(AGYW, FSW, MEN, CFSW, PREG)[2],
    MEN_rank = rank_populations(AGYW, FSW, MEN, CFSW, PREG)[3],
    CFSW_rank = rank_populations(AGYW, FSW, MEN, CFSW, PREG)[4],
    PREG_rank = rank_populations(AGYW, FSW, MEN, CFSW, PREG)[5]
  ) %>%
  ungroup()

# Join with parameter values
df_diagnosis <- impact_ranked |>
  filter(test_strat == "Diagnostic testing") |>
  left_join(post_samples, by = "lhs_num") 

df_screening <- impact_ranked |>
  filter(test_strat == "Screening") |>
  left_join(post_samples, by = "lhs_num") 

# ANALYSE RANKS ----

## Diagnosis
df_diagnosis |>
  select(lhs_num, ends_with("_rank")) |>
  pivot_longer(cols = ends_with("_rank"), 
               names_to = "population", 
               values_to = "rank") |>
  group_by(population, rank) |>
  summarise(count = n()) |>
  mutate(freq = count/sum(count)) |>
  arrange(population, rank) |> print(n=100)

df_diagnosis |>
  select(lhs_num, ends_with("_rank")) |>
  pivot_longer(cols = ends_with("_rank"), 
               names_to = "population", 
               values_to = "rank") |>
  group_by(population, rank) |>
  summarise(count = n()) |>
  mutate(freq = count/sum(count)) |>
  arrange(population, rank) |> 
  filter(rank == 1)

df_diagnosis |>
  filter(PREG_rank == 1)
  

# proportion ranked according to median
df_diagnosis |>
  select(lhs_num, ends_with("_rank")) |>
  mutate(original = case_when(CFSW_rank == 1 & FSW_rank == 2 & MEN_rank == 3 & PREG_rank == 4 & AGYW_rank == 5 ~ "yes",
                              TRUE ~ "no")) |>
  group_by(original) |>
  summarise(n= n())

# plots ranks
df_diagnosis |>
  select(lhs_num, ends_with("_rank")) |>
  pivot_longer(cols = ends_with("_rank"), 
               names_to = "population", 
               values_to = "rank") |>
  ggplot(aes(x = rank, fill = population)) +
  geom_bar() +
  facet_grid(~population)

df_diagnosis |> filter(FSW_rank == 3)

## Screening
df_screening |>
  select(lhs_num, ends_with("_rank")) |>
  pivot_longer(cols = ends_with("_rank"), 
               names_to = "population", 
               values_to = "rank") |>
  group_by(population, rank) |>
  summarise(count = n()) |>
  mutate(freq = count/sum(count)) |>
  arrange(population, rank)

df_screening |>
  select(lhs_num, ends_with("_rank")) |>
  pivot_longer(cols = ends_with("_rank"), 
               names_to = "population", 
               values_to = "rank") |>
  group_by(population, rank) |>
  summarise(count = n()) |>
  mutate(freq = count/sum(count)) |>
  arrange(population, rank) |> 
  filter(rank == 1)

# proportion ranked according to median
df_screening |>
  select(lhs_num, ends_with("_rank")) |>
  mutate(original = case_when(FSW_rank == 1 & CFSW_rank == 2 & AGYW_rank == 3 & MEN_rank == 4 & PREG_rank == 5 ~ "yes",
                              TRUE ~ "no")) |>
  group_by(original) |>
  summarise(n= n())

# plot ranks
df_screening |>
  select(lhs_num, ends_with("_rank")) |>
  pivot_longer(cols = ends_with("_rank"), 
               names_to = "population", 
               values_to = "rank") |>
  ggplot(aes(x = rank, fill = population)) +
  geom_bar() +
  facet_wrap(~population)

# PRCC FOR RANK ----

# Get parameter names 
param_list <- names(post_samples)[1:(which(names(post_samples) == "lhs_num") - 1)]

# Function to analyze rank reversals using PRCC
rank_prcc <- function(data, params, col, n_bootstrap = 100) {
  
  # Create matrix X with parameter columns
  X <- as.data.frame(data[, params])
  
  # Rename columns to avoid problematic characters
  names(X) <- gsub("-", "_", params)
  
  # Create outcome vector - convert logical to numeric (0/1)
  y <- as.numeric(data[[col]])
  
  # Calculate PRCC
  prcc_results <- sensitivity::pcc(X = X, y = y, rank = TRUE, nboot = n_bootstrap)
  
  # Format results
  prcc_results$PRCC %>%
    rownames_to_column("parameter") %>%
    rename(prcc = original,
           prcc_lwr = `min. c.i.`,
           prcc_upr = `max. c.i.`) %>%
    mutate(
      sig = case_when(
        (prcc_lwr >= 0 & prcc_upr >= 0) | 
          (prcc_lwr <= 0 & prcc_upr <= 0) ~ "TRUE",
        TRUE ~ "FALSE"),
      pop = col
    )
}

# DIAGNOSIS
# Loop through all reversal columns
diagnosic_results <- list()
i <- 1

col <- c("AGYW_rank", "FSW_rank", "MEN_rank", "CFSW_rank", "PREG_rank")

for (r in col) {
  # Calculate PRCC for this scenario
  result <- rank_prcc(df_diagnosis, param_list, r)
  # Store result in list
  diagnosic_results[[i]] <- result
  i <- i + 1
}

# SCREENING
# Loop through all reversal columns
screening_results <- list()
i <- 1

for (r in col) {
  result <- rank_prcc(df_screening, param_list, r)
  screening_results[[i]] <- result
  i <- i + 1
}

# Look at significant parameters for each reversal type
bind_rows(diagnosic_results) |>
  filter(sig == "TRUE") |>
  select(parameter, prcc, pop) |>
  pivot_wider(names_from = pop, values_from = prcc) |>
  print(n=100)

bind_rows(screening_results) |>
  filter(sig == "TRUE") |>
  select(parameter, prcc, pop) |>
  pivot_wider(names_from = pop, values_from = prcc) |>
  print(n=100)

bind_rows(screening_results) |>
  filter(sig == "TRUE") |>
  select(parameter, prcc, prcc_lwr, prcc_upr, pop) |>
  filter(pop == "FSW_rank") |>
  arrange(prcc)

bind_rows(screening_results) |>
  filter(sig == "TRUE") |>
  select(parameter, prcc, prcc_lwr, prcc_upr, pop) |>
  filter(pop == "CFSW_rank") |>
  arrange(prcc)

# PLOTS ----

plot_diag <- 
  bind_rows(diagnosic_results) |>
  left_join(lhs_input |> select(var_strat, var_ggplot) |> mutate(var_strat = gsub("-", "_", var_strat)), 
            by = c("parameter"="var_strat")) |>
  mutate(parameter = factor(parameter, levels = gsub("-", "_", order)),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])),
         pop = str_replace(pop, "_rank", ""),
         pop = factor(pop, levels = c("CFSW","FSW", "MEN","PREG","AGYW"),
                      labels = c("CFSW", "FSW", "Men", "Pregnant", "AGYW"))) |>
  # categorise as negligible, weak, moderate, strong 
  mutate(cat = case_when(abs(prcc) >= 0.6 ~ "Strong",
                         abs(prcc) >= 0.4 ~ "Moderate",
                         abs(prcc) >= 0.2 ~ "Weak",
                         TRUE ~ "Negligible"),
         cat = factor(cat, levels = c("Negligible", "Weak", "Moderate", "Strong")))

plot_screen <- 
  bind_rows(screening_results) |>
  left_join(lhs_input |> select(var_strat, var_ggplot) |> mutate(var_strat = gsub("-", "_", var_strat)), 
            by = c("parameter"="var_strat")) |>
  mutate(parameter = factor(parameter, levels = gsub("-", "_", order)),
         var_ggplot = factor(var_ggplot, levels = unique(lhs_input$var_ggplot[match(order, lhs_input$var_strat)])),
         pop = str_replace(pop, "_rank", ""),
         pop = factor(pop, levels = c("CFSW","FSW", "MEN","PREG","AGYW"),
                      labels = c("CFSW", "FSW", "Men", "Pregnant", "AGYW"))) |>
  # categorise as negligible, weak, moderate, strong 
  mutate(cat = case_when(abs(prcc) >= 0.6 ~ "Strong",
                         abs(prcc) >= 0.4 ~ "Moderate",
                         abs(prcc) >= 0.2 ~ "Weak",
                         TRUE ~ "Negligible"),
         cat = factor(cat, levels = c("Negligible", "Weak", "Moderate", "Strong")))

colours <- c("grey70", "#B8E5CC", "#48B87C", "#2A8B65")

p1 <- plot_diag |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6")) +
  labs(x = "", y = "", tag = "A", fill = "")


p2 <- plot_screen |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6"),
                    guide = "none") +
  labs(x = "", y = "", tag = "B", fill = "")

p1 / p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(unit(c(t=-10,r=0,b=5,l=-5), "cm")),
        plot.margin = margin(unit(c(t=5,r=2,l=0,b=0),"cm")))

ggsave(".././plots/prcc_ranks.png", width = 18, height = 25, unit = "cm", dpi = 700)  

plot_diag |>
  filter(parameter == "zeta_1") |>
  select(pop, prcc, prcc_lwr, prcc_upr)

# PLOT SEPARATE ----

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

p1 <- plot_diag |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6")) +
  labs(x = "", y = "", tag = "A", fill = "")


p2 <- plot_screen |>
  ggplot(aes(y = forcats::fct_rev(var_ggplot))) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_vline(xintercept = -0.5, linetype = "dotted", linewidth = 0.3, colour = "grey30") + 
  geom_bar(aes(x = prcc, fill = cat), stat = "identity", width = 0.7) +
  geom_linerange(aes(xmin = prcc_lwr, xmax = prcc_upr), linewidth = 0.2) +
  facet_grid(~pop) +
  mytheme +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0.05)) +  
  scale_y_discrete(labels = function(x) parse(text=x)) +
  scale_fill_manual(values = colours, breaks = c("Negligible", "Weak", "Moderate", "Strong"),
                    labels = c("Very weak: |PRCC| < 0.2", "Weak: 0.2 ≤ |PRCC| < 0.4", 
                               "Moderate: 0.4 ≤ |PRCC| < 0.6", "Strong: |PRCC| ≥ 0.6")) +
  labs(x = "", y = "", tag = "B", fill = "")

ggsave(".././plots/prcc_ranks_diag.png", p1, width = 18, height = 18, unit = "cm", dpi = 700)  
ggsave(".././plots/prcc_ranks_screen.png", p2, width = 18, height = 18, unit = "cm", dpi = 700)  
