##### Appendix 1: Natural regeneration preliminary study ####
## Summary statistics and basic visualizations: ##
## ## ###

## Library ####
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(lme4)
library(DHARMa)
library(glmmTMB)  
library(broom.mixed)


# Read raw data
nat_regen <- read.csv(here("Raw_data","Nat_regen.csv"))
EIV_dat<- read.csv(here("Raw_data","NR_EIV.csv"))

# Step 1: Calculate total and max germination per row
NR_summary <- nat_regen %>%
  # Convert from wide to long format to easily summarize across seasons
  pivot_longer(
    cols = Autumn_2022:Spring_2025,
    names_to = "Season",
    values_to = "Germination_Count"
  ) %>%
  # Group by each unique plot to calculate summary stats
  group_by(Site, Plot, Mother_tree) %>%
  summarise(
    total_germ = sum(Germination_Count, na.rm = TRUE),
    max_germ = max(Germination_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Calculate germination rates using the defined parameter
  mutate(
    total_germ_rate = total_germ / 50,
    max_germ_rate = max_germ / 50
  )

# Step 2: Average per maternal tree (across all their plots)
NR_maternal_summary <- NR_summary %>%
  group_by(Site, Mother_tree) %>%
  summarise(
    mean_total_germ_rate = mean(total_germ_rate, na.rm = TRUE),
    mean_max_germ_rate = mean(max_germ_rate, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Summary by site — mean and standard error
NR_maternal_summary %>%
  group_by(Site) %>%
  summarise(
    mean_total_germ = mean(mean_total_germ_rate),
    se_total_germ = sd(mean_total_germ_rate) / sqrt(n()),
    mean_max_germ = mean(mean_max_germ_rate),
    se_max_germ = sd(mean_max_germ_rate) / sqrt(n()),
    .groups = "drop"
  )

# Step 4: Plotting mean_max_germ_rate by maternal tree and site
NR_p<-ggplot(NR_maternal_summary, aes(x = Mother_tree, y = mean_max_germ_rate, fill = Site)) +
  geom_col(position = position_dodge(width = 0.8)) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    aes(group = Site),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  labs(
    title = "Mean Maximum Germination Rate by Maternal Tree and Site",
    x = "Maternal Tree",
    y = "Germination Rate"
  ) +
  theme_minimal() 

ggsave(
  filename = here("Results", "NR_germination_overall.png"),
  plot = NR_p,
  width = 8, height = 6, dpi = 300,  bg = "white" )
## test for differences in germination counts for maternal trees####

output_file <- here("Results", "NR_maternal_germ_model.txt")
sink(output_file)

cat("### Maternal Tree Germination Variation LRT Results\n\n")

for (site_i in unique(NR_summary$Site)) {
  cat("\n--- Site:", site_i, "---\n\n")
  
  site_data <- filter(NR_summary, Site == site_i)
  
  # Full model with Plot/Mother_tree random effects
  model_full <- glmmTMB(
    cbind(total_germ, 50 - total_germ) ~ (1 | Plot/Mother_tree),
    family = binomial(link = "logit"),
    data = site_data
  )
  
  # Reduced model without maternal tree random effect (only Plot)
  model_reduced <- glmmTMB(
    cbind(total_germ, 50 - total_germ) ~ (1 | Plot),
    family = binomial(link = "logit"),
    data = site_data
  )
  
  cat("Full Model Summary:\n")
  print(summary(model_full))
  
  cat("\nReduced Model Summary:\n")
  print(summary(model_reduced))
  
  cat("\nVariance Components (Full Model):\n")
  print(VarCorr(model_full))
  
  # Likelihood ratio test
  lrt <- anova(model_reduced, model_full)
  cat("\nLikelihood Ratio Test comparing Reduced vs Full Model:\n")
  print(lrt)
  
  ## Residual Diagnostics
  cat("\n--- Model Diagnostics (DHARMa Simulated Residuals) ---\n")
  sim_res <- simulateResiduals(fittedModel = model_full, n = 1000)
  
  ## Test overdispersion
  cat("\nOverdispersion Test:\n")
  print(testOverdispersion(sim_res))
  
  ## Test zero-inflation
  cat("\nZero-Inflation Test:\n")
  print(testZeroInflation(sim_res))
  
  ## Uniformity test (residuals)
  cat("\nUniformity Test:\n")
  print(testUniformity(sim_res))
  
  ## Dispersion and residual outliers
  cat("\nOutlier Test:\n")
  print(testOutliers(sim_res))
  cat("\n------------------------------------------------------------\n")
}

sink()



#### EIV values : ####
## Species richness per plot: 
EIV_dat %>% 
  mutate(species_id = paste(Genus, Species)) %>%
  group_by(Site, Plot) %>%
  summarise(
    species_richness = n_distinct(species_id),
    .groups = "drop"
  )

## take means per plot : 
EIV_means <- EIV_dat %>%
  group_by(Site, Plot) %>%
  summarise(
    mean_L = mean(L, na.rm = TRUE),
    mean_F = mean(F, na.rm = TRUE),
    mean_R = mean(R, na.rm = TRUE),
    mean_N = mean(N, na.rm = TRUE),
    mean_S = mean(S, na.rm = TRUE),
    .groups = "drop"
  )

## test for differences in mean values across plots: 
ellenberg_traits <- c(
  "mean_L" = "Light",
  "mean_F" = "Moisture",
  "mean_R" = "pH",
  "mean_N" = "Nutrients",
  "mean_S" = "Salinity"
)

output_file <- here("Results", "All_EIV_models.txt")

# Send all printed output to this file
sink(output_file)

# Loop through sites and Ellenberg traits
for (site in unique(EIV_means$Site)) {
  cat("\n==============================\n")
  cat("Site:", site, "\n")
  site_data <- filter(EIV_means, Site == site)
  
  for (trait in names(ellenberg_traits)) {
    label <- ellenberg_traits[[trait]]
    
    formula <- as.formula(paste(trait, "~ Plot"))
    model <- aov(formula, data = site_data)
    
    cat("\n--- Trait:", label, "(", trait, ") ---\n")
    cat("ANOVA results:\n")
    print(summary(model))
    
    # Shapiro-Wilk test of residuals (normality)
    shapiro_p <- shapiro.test(resid(model))$p.value
    cat("Shapiro-Wilk normality p =", round(shapiro_p, 4), "\n")
    
    # Residual summary
    res <- resid(model)
    cat("Residuals: min =", round(min(res), 3), 
        "| max =", round(max(res), 3), 
        "| mean =", round(mean(res), 3), "\n")
  }
}

sink()

### plot EIV values
EIV_long <- EIV_dat %>%
  pivot_longer(
    cols = c(L, F, R, N, S),
    names_to = "EIV",
    values_to = "Value"
  ) %>%
  mutate(EIV = dplyr::recode(EIV,
                      "L" = "Light",
                      "F" = "Moisture",
                      "R" = "pH",
                      "N" = "Nutrients",
                      "S" = "Salinity"))



ggplot(EIV_long, aes(x = factor(Plot), y = Value, fill = factor(Plot))) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  facet_wrap(~ Site + EIV, scales = "free_y", ncol = 5) +  # Facet by Site and EIV, vertical layout
  labs(
    title = "Ellenberg Indicator Values by Plot within Sites",
    x = "Plot",
    y = "Ellenberg Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Hide legend since fill = Plot is obvious on x-axis
  )


### Test for EIV predicting germ: ####
## join datasets: 
data_model <- NR_summary %>%
  left_join(EIV_means, by = c("Site", "Plot"))

# Fit binomial GLMM with nesting and EIV traits
model <- glmer(
  cbind(total_germ, 50 - total_germ) ~ Site + mean_L + mean_F + mean_R + mean_N + mean_S +
    (1 | Site/Plot/Mother_tree),
  family = binomial(link = "logit"),
  data = data_model,
  control = glmerControl(optimizer = "bobyqa")
)

# Output summary and diagnostics to file
sink(here("Results","NR_germ_EIV_model.txt"))
cat("### Germination Model Summary ###\n\n")
print(summary(model))

cat("\n\n### Variance Components ###\n\n")
print(VarCorr(model), comp="Std.Dev.")

cat("\n\n### Model Diagnostics (DHARMa) ###\n\n")
sim_res <- simulateResiduals(model)
plot(sim_res)
testDispersion(sim_res)
testZeroInflation(sim_res)
testUniformity(sim_res)
testOutliers(sim_res)
sink()