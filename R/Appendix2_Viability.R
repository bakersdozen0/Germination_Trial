#### Library ####
## Please install the following required packages if not already installed in
## your environment. The library is shared across all scritps for consistency:
## Reproducibility: 
## NB: here() commands depend on Rproject being saved in same directory as data
library(here)
## Data management and visualization tools:
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lubridate)
library(broom)
## For zero-inflated modelling and model comparison/validation:
library(lme4)
library(glmmTMB)
library(DHARMa)
## For time-to event modelling:
library(drcte)
library(drcSeedGerm)


#### Appendix 2: Viability analyses ####
#### Create or Load full_fam_germ

full_fam_germ_path <- here("Proc_data", "full_fam_germ.csv")

if (!file.exists(full_fam_germ_path)) {
  cat("Processed full_fam_germ.csv not found. Generating...\n")
  
  fam_germ <- plate_summary %>%
    group_by(Fam) %>%
    summarise(
      total_germinated = sum(total_germinated, na.rm = TRUE),
      total_seeds = sum(Starting_number, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(fam_fam_germ_percent = total_germinated / total_seeds)
  
  full_fam_germ <- fam_germ %>%
    left_join(master, by = "Fam") %>%
    left_join(TZV, by = "Fam") %>%
    left_join(XRV, by = "Fam")
  
  write.csv(full_fam_germ, full_fam_germ_path, row.names = FALSE)
  
} else {
  
  cat("Found full_fam_germ.csv. Loading from file...\n")
  full_fam_germ <- read.csv(full_fam_germ_path)
}


#### estimate numbers of seeds/ grams of berry collected: 
seed_density_summary <- full_fam_germ %>%
  mutate(
    # Estimate number of filled seeds in the batch
    est_filled_seeds = Weight_of_filled_batch / Average_filled_seed_weight,
    
    # Estimate number of empty seeds in the batch
    est_empty_seeds = Weight_of_empty_batch / Average_empty_seed_weight,
    
    # Total estimated seeds in batch
    est_total_seeds = est_filled_seeds + est_empty_seeds,
    
    # Seeds per gram of total batch weight
    seeds_per_gram = est_total_seeds / Batch_weight
  ) %>%
  summarise(
    avg_filled_seeds_per_gram = mean(est_filled_seeds / Batch_weight, na.rm = TRUE),
    avg_empty_seeds_per_gram = mean(est_empty_seeds / Batch_weight, na.rm = TRUE),
    avg_total_seeds_per_gram = mean(seeds_per_gram, na.rm = TRUE)
  )

print(seed_density_summary)

full_fam_germ<-full_fam_germ %>% 
  mutate(filled_prop= (Average_filled_seed_weight*Weight_of_filled_batch)/((Average_filled_seed_weight*Weight_of_filled_batch)+(Average_empty_seed_weight*Weight_of_empty_batch)) )

## Spearman's correlations among germination and viability, and "like" viability meas.  ####
## first simply run correlations - how are the viability methods related to one 
## another and to actual overall family germination? 

full_fam_germ<-full_fam_germ %>% 
  mutate(across(c(Parasite_load_percent, FCT_percent,FCT_green_berries), ~ . / 100))

# Scale proportions to (0,1)
scale_proportion <- function(x) {
  scaled <- (x / 100) * 0.99 + 0.005
  return(scaled)
}

full_fam_germ <- full_fam_germ %>%
  mutate(
    FCT_percent_prop = scale_proportion(FCT_percent),
    FCT_green_berries_prop = scale_proportion(FCT_green_berries),
    log_FE_ratio = log(Filled_fraction)
  )

# Remove outliers from F:E ratio
fe_iqr <- quantile(full_fam_germ$Filled_fraction, probs = c(0.25, 0.75), na.rm = TRUE)
iqr <- diff(fe_iqr)
fe_bounds <- fe_iqr + c(-1.5, 1.5) * iqr

full_fam_germ <- full_fam_germ %>%
  filter(Filled_fraction >= fe_bounds[1], Filled_fraction <= fe_bounds[2])


### 
sink(here("Results","Viability_correlations.txt"))

cat( "==== All viability correlations ====")
cat( "= Firstly, correlations with germination: =")
cor.test(full_fam_germ$germ_percent, full_fam_germ$FCT_percent, method = "spearman")
cor.test(full_fam_germ$germ_percent, full_fam_germ$FCT_green_berries, method = "spearman")
cor.test(full_fam_germ$germ_percent, full_fam_germ$Filled_fraction, method = "spearman")
cor.test(full_fam_germ$germ_percent, full_fam_germ$TZ_Viability, method = "spearman")
cor.test(full_fam_germ$germ_percent, full_fam_germ$XR_Viability, method = "spearman")

# Correlations between filled vs. empty seed measures:
cat( "= Secondly, correlations between measures of filled v. empty seeds: =")
cor.test(full_fam_germ$FCT_percent, full_fam_germ$Filled_fraction, method = "spearman")
cor.test(full_fam_germ$FCT_percent, full_fam_germ$FCT_green_berries_prop, method = "spearman")
cor.test(full_fam_germ$FCT_green_berries_prop, full_fam_germ$Filled_fraction, method = "spearman")

cat( "= Lastly, correlations betweent measures of metabolic viability: =")
#Correlations between viable and nonviable seeds
cor.test(full_fam_germ$TZ_Viability, full_fam_germ$XR_Viability, method = "spearman")

sink()

### All viability models ####
# Prepare data 


# Models :: 
# Reported in paper are these specifications for each seed lot quality/ viability metric: 
# success/failure (beta-family) OR proportion / ratio (lm family) ~ Pop, OR  ~ Region +(1|Region/Pop)

## NB: Reference currently BT, reset to FC (closest value to overall average and median)
full_fam_germ <- full_fam_germ %>%
  mutate(across(c(Region, Pop, Fam), as.factor)) %>%
  mutate(Pop = relevel(Pop, ref = "FC"))


# FCT viability model
mod_fct <- glmmTMB(FCT_percent_prop ~ Pop,
                   family = beta_family(),
                   data = full_fam_germ)

# Green berries (proportion) model
mod_green <- glmmTMB(FCT_green_berries_prop ~ Pop,
                     family = beta_family(),
                     data = full_fam_germ)

# Filled:Empty Ratio (non-transformed) model
mod_fe <- lm(Filled_fraction ~ Pop, 
             data = full_fam_germ)

# Filled:Empty Ratio (log-transformed) model
mod_logfe <- lm(log_FE_ratio ~ Pop, 
                data = full_fam_germ)

## TZ by region and Pop, and XR by pop have signifcant DHARMA residuals as binomial families, change to beta-binomial to account for 
## overdispersion: 

# TZ viability (count data)
mod_TZ <- glmmTMB(cbind(Stained, Sum_tested - Stained) ~ Pop,
                  family = betabinomial(link="logit"), data = full_fam_germ)

# XR viability (count data)
mod_XR <- glmmTMB(cbind(Viable, X_ray_sample_size - Viable) ~ Pop,
                  family =  betabinomial(link="logit"), data = full_fam_germ)

# --- Region-level Viability Models ---

# FCT (Proportion) by Region
mod_fct_region <- glmmTMB(FCT_percent_prop ~ Region + (1|Region:Pop),
                          family = beta_family(),
                          data = full_fam_germ)

# XR Viability (Counts) by Region
mod_XR_region <- glmmTMB(cbind(Viable, X_ray_sample_size - Viable) ~ Region + (1|Region:Pop),
                         family = betabinomial(link="logit"), 
                         data = full_fam_germ)

# TZ Viability (Counts) by Region
mod_TZ_region <- glmmTMB(cbind(Stained, Sum_tested - Stained) ~ Region + (1|Region:Pop),
                         family = betabinomial(link="logit"), 
                         data = full_fam_germ)

## Write list of models to iterate over for saving results: 
models <- list(
  mod_fct = mod_fct,
  mod_green = mod_green,
  mod_fe = mod_fe,
  mod_logfe=mod_logfe,
  mod_TZ = mod_TZ,
  mod_XR = mod_XR,
  mod_fct_region=mod_fct_region,
  mod_XR_region=mod_XR_region,
  mod_TZ_region
)

### Function for saving model results using sink()
for (name in names(models)) {
  model <- models[[name]]
  
  cat("Running diagnostics and saving outputs for:", name, "\n")
  
  # Open text file sink for summary + diagnostics
  sink(here("Results", paste0(name, "_by_Pop_summary.txt")))
  
  cat("===== Model Summary: ", name, " =====\n\n")
  print(summary(model))
  
  # Only run DHARMa diagnostics for glmmTMB or similar models
  if (inherits(model, "glmmTMB") || inherits(model, "glmerMod")) {
    cat("\n\n===== DHARMa Diagnostics =====\n")
    simres <- simulateResiduals(model)
    print(testDispersion(simres))
    print(testUniformity(simres))
    print(testZeroInflation(simres))
  } else {
    cat("\n\n(No DHARMa diagnostics available for this model type)\n")
  }
  
  sink()
  
  # Save diagnostic plot if applicable
  if (exists("simres")) {
    png(here("Results", paste0(name, "_residuals.png")), width = 800, height = 600)
    plot(simres)
    dev.off()
    rm(simres)  # remove simres to avoid confusion in next iteration
  }
}
