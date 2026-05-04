#### R code for simple transformations and analyses, following those done in
#### Baker et. al (2025, in prep) : see section headers for details  
#### created 19/6/25 by JB, in preparation of submitting data to EIDC, and paper
#### to publications. Includes transformations to long format, standardizing 
#### series dates (see paper), calculating descriptive statistics, both zero
#### inflated modelling and time-to-event modelling, and simple visualizations
#### Please see metadata in EIDC deposit for a description of the data types
#### and general analyses. 


#### Library ####
## Please install the following required packages if not already installed in
## your environment.
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

#### Data import and transformations: ####
## Please see below for code to conditionally load or created processed files:
## NB: these commands assume that the names/ file structure of the files from the 
## EIDC submission have not been changed.

## Raw data import: 
## Read in .csvs for collection data, plate checks over time:
master<-read.csv(here("Raw_data","Germination_master.csv"))
PCs<-read.csv(here("Raw_data","Plate_checks.csv"))
## X-ray and TZ viability testing: 
TZV<-read.csv(here("Raw_data","TZ_Viability.csv"))
XRV<-read.csv(here("Raw_data","XR_Viability.csv"))

str(PCs)
## Processed data: 
## NB: These processed data were not included in EIDC submission, 
## use these commands to import them after they are generated, following the script below
#PCs_long <- read.csv(here("Proc_data", "PCs_long.csv"))%>% mutate(across(c(Region, Pop, Fam, Treatment, Pretreatment), as.factor))
#plate_summary <- read.csv(here("Proc_data", "plate_summary.csv"))%>% mutate(across(c(Region, Pop, Fam, Treatment, Pretreatment), as.factor))
#pop_DS <- read.csv(here("Results", "Treat_DS.csv"))
#tte_data <- read.csv(here("Proc_data", "tte_data.csv"))%>% mutate(across(c(Region, Pop, Fam, Treatment, Pretreatment), as.factor))

## see appendix 2 for code to generate full_fam_germ
#full_fam_germ <- read.csv(here("Proc_data", "full_fam_germ.csv"))%>% mutate(across(c(Region, Pop, Fam, Treatment, Pretreatment), as.factor))


## These typically read into R ok, we will need to adjust data types for analyses
## (see below), but one general fix we do need to make is how R reads the time-series data
## Please see metadata for an explanation of these data 

#### STEP 1: Create or Load Processed data for all analyses: ####

# Define the path for the processed file
processed_pcs_long_path <- here("Proc_data", "PCs_long.csv")

# Check if the processed file exists. If not, create it.
## NB: this will throw non-lethal warnings, this is just about introducing NA's
if (!file.exists(processed_pcs_long_path)) {
  
  cat("Processed file not found. Generating PCs_long.csv from raw data...\n")
  #First remove "X"'s from date columns:
  colnames(PCs) <- gsub("^X", "", colnames(PCs))
  #Convert check data to numeric:
  PCs <- PCs %>%
    mutate(across(12:49, as.numeric))
  #convert to long format and perform transformations
  PCs_long <- PCs %>%
    pivot_longer(cols = 12:49,
                 names_to = "Check_Date",
                 values_to = "Germinated") %>%
    mutate(
      # Change date columns to date format
      Check_Date = as.Date(Check_Date, format = "%d.%m.%y"),
      Stratification_start_date = as.Date(Stratification_start_date, format = "%d.%m.%y"),
      # Convert key columns to factors
      across(c(Region, Pop, Fam, Treatment, Pretreatment), as.factor)
    ) %>%
    # Add a data validation/sanity check
    filter(!is.na(Germinated) & Germinated <= Starting_number)
  # Write the newly created file for next time
  write.csv(PCs_long, processed_pcs_long_path, row.names = FALSE)
} else {
  cat("Found processed_pcs_long.csv. Loading from file...\n")
  PCs_long <- read.csv(processed_pcs_long_path)
  # Ensure correct data types after loading from CSV
  PCs_long <- PCs_long %>%
    mutate(
      Check_Date = as.Date(Check_Date),
      Stratification_start_date = as.Date(Stratification_start_date),
      across(c(Region, Pop, Fam, Treatment, Pretreatment), as.factor)
    )
}

## Check: Check dates should span Nov. 2023- May 2025
str(PCs_long)
unique(PCs_long$Check_Date)

#### STEP 2: Transformations on the PCs_long Object
## Standardize dates by start of each treatment (account for Wave);
## NB: The terminology "wave" is used here, "series" is used in the paper. 
## these terms are synonymous. The trial was started over 5 "waves"/"series"
PCs_long <- PCs_long %>% 
  mutate(Days_since_start = as.numeric(Check_Date - Stratification_start_date)) %>% 
  mutate(standardized_weeks_start = Days_since_start / 7)

#Standardize by end dates of each treatment: (account for wave and treatment)
treatment_durations <- data.frame(
  Treatment = c("Control", "W", "W->C", "C"),
  strat_duration = c(0, 8, 23, 15)
)

#add length of treatments
PCs_long <- merge(PCs_long, treatment_durations, by = "Treatment", all.x = TRUE)

#add end of treatments as a date to calculate number of weeks
PCs_long$end_stratification <- PCs_long$Stratification_start_date + (PCs_long$strat_duration * 7)

#Calculate the number of weeks from the end of stratification to the observation date
PCs_long$standardized_weeks_end <- as.numeric(difftime(PCs_long$Check_Date,
                                                       PCs_long$end_stratification, 
                                                       units = "weeks"))
# Check the final dataset
head(PCs_long)
str(PCs_long)
#### Calculate plate_summary (seed and germination counts) per replicate) :
### join with family-level x-ray and TZ data for later 

plate_summary_path <- here("Proc_data", "plate_summary.csv")

### 
if (!file.exists(plate_summary_path)) {
  
  cat("Processed plate_summary.csv not found. Generating from PCs_long...\n")
  
  plate_summary <- PCs_long %>%
    filter(!is.na(Germinated)) %>%
    group_by(Plate_ID) %>%
    summarise(
      total_germinated = sum(Germinated, na.rm = TRUE),
      Starting_number = first(Starting_number),
      Region = first(Region),
      Pop = first(Pop),
      Fam = first(Fam),
      Treatment = first(Treatment),
      Pretreatment = first(Pretreatment),
      .groups = "drop"
    ) %>%
    left_join(XRV %>% dplyr::select(Fam, XR_Viability), by = "Fam") %>%
    mutate(
      Germ_Percent = total_germinated / Starting_number,
      No_Germination = total_germinated == 0
    )
  write.csv(plate_summary, plate_summary_path, row.names = FALSE)
  
} else {
  
  cat("Found plate_summary.csv. Loading from file...\n")
  plate_summary <- read.csv(plate_summary_path)
  
  # Ensure proper types (if needed)
  plate_summary <- plate_summary %>%
    mutate(across(c(Region, Pop, Fam, Treatment, Pretreatment), as.factor))
}
str(plate_summary)

#### STEP 4:Create or Load Pop or Treat Summary (Pop_DS.csv/ Treat_DS.csv)
## edit here(), if() , pop_DS and group_by() to do this by treatment/Pop/Region

## define function to calculate DS: 
## grouping_variable either "Pop" or "Treatment" (or optionally Pretreatment, Fam etc.)
calculate_DS <- function(data, grouping_variable) {
  data %>%
    group_by({{ grouping_variable }}) %>%
    summarise(
      Mean_Germ_Percent = mean(Germ_Percent, na.rm = TRUE),
      Median_Germ_Percent = median(Germ_Percent, na.rm = TRUE),
      SD_Germ_Percent = sd(Germ_Percent, na.rm = TRUE),
      N_Plates = n(),
      SE_Germ_Percent = SD_Germ_Percent / sqrt(N_Plates),
      CV_Germ_Percent = SD_Germ_Percent / Mean_Germ_Percent,
      Prop_Zero_Germ = mean(No_Germination, na.rm = TRUE),
      .groups = "drop"
    )
}


## apply function: 
## NB: X-ray viability only applied to Pop-level summary: 
pop_ds_path <- here("Results", "Pop_DS.csv")

if (!file.exists(pop_ds_path)) {
  cat("Generating Pop_DS.csv from plate_summary...\n")
  
  # Step 1: Calculate core germination stats using the function
  pop_DS <- calculate_DS(plate_summary, Pop)
  
  # Step 2: Add Population-specific data (mean X-ray viability)
  # Create a mapping from family to population
  fam_to_pop <- PCs_long %>%
    dplyr::select(Fam, Pop) %>%
    distinct()
  
  # Calculate mean viability per population
  xr_viability_pop <- XRV %>%
    left_join(fam_to_pop, by = "Fam") %>%
    group_by(Pop) %>%
    summarise(Mean_XR_Viability = mean(XR_Viability, na.rm = TRUE), .groups = "drop")
  
  # Join the viability data to the main summary table
  pop_DS <- pop_DS %>%
    left_join(xr_viability_pop, by = "Pop")
  
  # Step 3: Write the final file
  write.csv(pop_DS, pop_ds_path, row.names = FALSE)
  
} else {
  cat("Found Pop_DS.csv. Loading from file...\n")
  pop_DS <- read.csv(pop_ds_path)
}

treatment_ds_path <- here("Results", "Treatment_DS.csv")

if (!file.exists(treatment_ds_path)) {
  cat("Generating Treatment_DS.csv from plate_summary...\n")
  
  # Step 1: Calculate stats using the function and save to a new object
  treatment_DS <- calculate_DS(plate_summary, Treatment)
  
  # Step 2: Write the file (no extra data to join for this grouping)
  write.csv(treatment_DS, treatment_ds_path, row.names = FALSE)
  
} else {
  cat("Found Treatment_DS.csv. Loading from file...\n")
  treatment_DS <- read.csv(treatment_ds_path)
}

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
    mutate(germ_percent = total_germinated / total_seeds)
  
  full_fam_germ <- fam_germ %>%
    left_join(master, by = "Fam") %>%
    left_join(TZV, by = "Fam") %>%
    left_join(XRV, by = "Fam")
  
  write.csv(full_fam_germ, full_fam_germ_path, row.names = FALSE)
  
} else {
  
  cat("Found full_fam_germ.csv. Loading from file...\n")
  full_fam_germ <- read.csv(full_fam_germ_path)
}


# Add columns to estimate number of seeds needed to get 100 germinated ##
# Step 1: Calculate filled fraction and seeds needed
fam_seed_requirements <- full_fam_germ %>%
  mutate(
    # Estimate number of filled and empty seeds using weight and average seed weight
    est_filled = ifelse(!is.na(Weight_of_filled_batch) & !is.na(Average_filled_seed_weight) & Average_filled_seed_weight > 0,
                        Weight_of_filled_batch / Average_filled_seed_weight,
                        NA_real_),
    
    est_empty = ifelse(!is.na(Weight_of_empty_batch) & !is.na(Average_empty_seed_weight) & Average_empty_seed_weight > 0,
                       Weight_of_empty_batch / Average_empty_seed_weight,
                       NA_real_),
    
    # Weight-based filled_fraction
    filled_fraction_weighted = ifelse(!is.na(est_filled) & !is.na(est_empty),
                                      est_filled / (est_filled + est_empty),
                                      NA_real_),
    
    # Estimate number of sown filled seeds to get 100 germinants
    sown_filled_seeds = ifelse(germ_percent > 0,
                               ceiling(100 / germ_percent),
                               NA_real_),
    
    # Estimate number of total seeds from full seedlot to sow (using weighted filled_fraction)
    seedlot_seeds_weighted = ifelse(germ_percent > 0 & !is.na(filled_fraction_weighted),
                                    ceiling(100 / (germ_percent * filled_fraction_weighted)),
                                    NA_real_)
  ) %>%
  dplyr::select(Fam, germ_percent, sown_filled_seeds,
                filled_fraction_weighted, seedlot_seeds_weighted)# Step 2: Summary statistics across families


seed_requirements_summary <- fam_seed_requirements %>%
  summarise(
    min_sown_filled_seeds = min(sown_filled_seeds, na.rm = TRUE),
    max_sown_filled_seeds = max(sown_filled_seeds, na.rm = TRUE),
    avg_sown_filled_seeds = mean(sown_filled_seeds, na.rm = TRUE),
    
    min_seedlot_seeds = min(seedlot_seeds_weighted, na.rm = TRUE),
    max_seedlot_seeds = max(seedlot_seeds_weighted, na.rm = TRUE),
    avg_seedlot_seeds = mean(seedlot_seeds_weighted, na.rm = TRUE)
  )

# View outputs
print(fam_seed_requirements)
View(seed_requirements_summary)

#### Summarize seed requirements by Population ##
# First, join the population information to the family-level data
fam_to_pop <- PCs_long %>%
  dplyr::select(Fam, Pop) %>%
  distinct()

fam_seed_requirements_with_pop <- fam_seed_requirements %>%
  left_join(fam_to_pop, by = "Fam")

# Then, group by Pop and calculate the average for each type of seed estimate
pop_seed_requirements_summary <- fam_seed_requirements_with_pop %>%
  group_by(Pop) %>%
  summarise(
    Avg_Sown_Filled_Seeds = mean(sown_filled_seeds, na.rm = TRUE),
    Avg_Seedlot_Seeds = mean(seedlot_seeds_weighted, na.rm = TRUE),
    .groups = "drop"
  )

# Print the data frame directly to the console
cat("Average estimated seeds needed per population:\n")
print(pop_seed_requirements_summary)

# Save the summary table to a CSV file in the Results folder
write.csv(pop_seed_requirements_summary, here("Results", "Pop_Seed_Requirements_Summary.csv"), row.names = FALSE)


# View outputs
print(fam_seed_requirements)
View(seed_requirements_summary)

##### for W only
# Step 1: Filter plate_summary for W treatment plates
plate_summary_W <- plate_summary %>%
  filter(Treatment == "W->C")

# Step 2: Calculate germ_percent per family within W treatment
fam_germ_W <- plate_summary_W %>%
  group_by(Fam) %>%
  summarise(
    germ_percent = mean(Germ_Percent, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Join with weight info (assuming full_fam_germ has weights and Fam)
fam_seed_requirements_W <- fam_germ_W %>%
  left_join(full_fam_germ %>% dplyr::select(Fam, Weight_of_filled_batch, Average_filled_seed_weight,
                                     Weight_of_empty_batch, Average_empty_seed_weight), by = "Fam") %>%
  mutate(
    est_filled = ifelse(!is.na(Weight_of_filled_batch) & !is.na(Average_filled_seed_weight) & Average_filled_seed_weight > 0,
                        Weight_of_filled_batch / Average_filled_seed_weight,
                        NA_real_),
    est_empty = ifelse(!is.na(Weight_of_empty_batch) & !is.na(Average_empty_seed_weight) & Average_empty_seed_weight > 0,
                       Weight_of_empty_batch / Average_empty_seed_weight,
                       NA_real_),
    filled_fraction_weighted = ifelse(!is.na(est_filled) & !is.na(est_empty),
                                      est_filled / (est_filled + est_empty),
                                      NA_real_),
    sown_filled_seeds = ifelse(germ_percent > 0,
                               ceiling(100 / germ_percent),
                               NA_real_),
    seedlot_seeds_weighted = ifelse(germ_percent > 0 & !is.na(filled_fraction_weighted),
                                    ceiling(100 / (germ_percent * filled_fraction_weighted)),
                                    NA_real_)
  ) %>%
  dplyr::select(Fam, germ_percent, sown_filled_seeds, filled_fraction_weighted, seedlot_seeds_weighted)

# Step 4: Summarize seed requirements across W treatment families
seed_requirements_summary_W <- fam_seed_requirements_W %>%
  summarise(
    min_sown_filled_seeds = min(sown_filled_seeds, na.rm = TRUE),
    max_sown_filled_seeds = max(sown_filled_seeds, na.rm = TRUE),
    avg_sown_filled_seeds = mean(sown_filled_seeds, na.rm = TRUE),
    
    min_seedlot_seeds = min(seedlot_seeds_weighted, na.rm = TRUE),
    max_seedlot_seeds = max(seedlot_seeds_weighted, na.rm = TRUE),
    avg_seedlot_seeds = mean(seedlot_seeds_weighted, na.rm = TRUE)
  )

# View summary
print(seed_requirements_summary_W)

#### TTE modelling: ####
## TTE modelling requires time series data for each plate, specially formatted for drtce
## These data first need to be standardized by the number of days since the stratification was 
## started. In other words, to account for different start dates and check dates Define (2 or) 3-week interval bins.
## What's more, drtce requires the same number of checks
## for every plate, regardless of treatment or wave (which is not the case, treatments are different lengths)
## so you need to populate a complete grid to represent data

## NB: I'll let the user decide what bins and cut offs to use. In the paper, 
## dates are binned by 3 week intervals (which is the frequency they were checked for most of the trial duration, see methods)
## and the cut off is 51 weeks, b/c the cumulative germination curve flattened by then
TTE_BIN_WIDTH_WEEKS <- 3 # The interval for grouping germination checks
TTE_MAX_WEEKS <- 51      # The maximum number of weeks to include in the model

tte_data_path <- here("Proc_data", "tte_data.csv")

if (!file.exists(tte_data_path)) {
  cat("Processed tte_data.csv not found. Generating from PCs_long...\n")
  
  max_observed_weeks <- ceiling(max(PCs_long$standardized_weeks_start, na.rm = TRUE))
  bin_edges <- seq(0, max_observed_weeks + TTE_BIN_WIDTH_WEEKS, by = TTE_BIN_WIDTH_WEEKS)
  bin_labels <- bin_edges[-length(bin_edges)]
  
  PCs_binned <- PCs_long %>%
    mutate(week_bin = cut(standardized_weeks_start, breaks = bin_edges, labels = bin_labels, right = FALSE, include.lowest = TRUE))
  
  plate_metadata <- PCs_binned %>%
    distinct(Plate_ID, ORDERALL, Region, Pop, Fam, Pretreatment, Treatment, Plate)
  
  bin_plate_grid <- expand_grid(Plate_ID = unique(PCs_binned$Plate_ID), week_bin = bin_labels) %>%
    mutate(week_bin = factor(week_bin, levels = bin_labels)) %>%
    left_join(plate_metadata, by = "Plate_ID")
  
  germ_summary <- PCs_binned %>%
    group_by(Plate_ID, week_bin) %>%
    summarise(total_seeds_germinated = sum(Germinated, na.rm = TRUE), .groups = "drop")
  
  complete_std_germ <- bin_plate_grid %>%
    left_join(germ_summary, by = c("Plate_ID", "week_bin")) %>%
    mutate(total_seeds_germinated = replace_na(total_seeds_germinated, 0))
  
  starting_seeds <- PCs %>%
    dplyr::select(Plate_ID, Starting_number)
  
  final_std_germ <- complete_std_germ %>%
    left_join(starting_seeds, by = "Plate_ID") %>%
    filter(as.numeric(as.character(week_bin)) <= TTE_MAX_WEEKS)
  
  tte_data_wide <- final_std_germ %>%
    pivot_wider(
      id_cols = c(Plate_ID, Region, Pop, Fam, Pretreatment, Treatment, Plate, ORDERALL, Starting_number),
      names_from = week_bin,
      values_from = total_seeds_germinated,
      values_fill = 0
    )
  
  # Identify the names of columns that are purely numeric (the week bins)
  count_column_names <- names(tte_data_wide)[which(sapply(names(tte_data_wide), function(x) grepl("^[0-9]+$", x)))]
  
  # Identify the names of all other columns to be used as metadata/ID
  id_column_names <- setdiff(names(tte_data_wide), count_column_names)
  
  montimes <- seq(0, TTE_MAX_WEEKS, by = TTE_BIN_WIDTH_WEEKS)
  
  tte_data <- melt_te(
    tte_data_wide,
    count_cols = all_of(count_column_names),
    treat_cols = all_of(id_column_names), # Explicitly define treatment/ID cols
    n.subjects = tte_data_wide$Starting_number,
    monitimes = montimes
  )
  
  tte_data <- tte_data %>%
    mutate(across(c(Treatment, Pretreatment, Region, Pop, Fam), as.factor))
  
  write.csv(tte_data, tte_data_path, row.names = FALSE)
} else {
  cat("Found tte_data.csv. Loading from file...\n")
  tte_data <- read.csv(tte_data_path) %>%
    mutate(across(c(Treatment, Pretreatment, Region, Pop, Fam), as.factor))
}


str(tte_data)

# Optional: sanity check on cumulative germination - should be 12, which is not an error (although maybe an accident)
# but a genuine and unique maximum value of germination in one plate: (due to "triplet" seed)
complete_std_germ %>%
  group_by(Plate_ID) %>%
  summarise(cumulative = sum(total_seeds_germinated, na.rm = TRUE)) %>%
  summarise(max_cumulative = max(cumulative)) %>%
  print()


## use loop to fit models for all data, as well as data by Region, Treatment and Region*Treatment:
fit_tte_model <- function(data, curveid = NULL, separate = FALSE, upperl = NULL) {
  drmte(
    count ~ timeBef + timeAf,
    fct = LL.3(),
    data = data,
    curveid = curveid,
    separate = separate,
    upperl = upperl
  )
}

## create list of modeling scenarios ## 
model_specs <- list(
  overall = list(curveid = NULL, separate = FALSE, upperl = NULL),
  treatment = list(curveid = quote(Treatment), separate = FALSE, upperl = c(NA, NA, NA, NA, 1, 1, 1, 1, NA, NA, NA, NA)),
  region = list(curveid = quote(Region), separate = FALSE, upperl = NULL),
  region_treatment = list(curveid = quote(Region:Treatment), separate = TRUE, upperl = NULL)
)

## apply function over list : NB: this will take several minutes
model_fits <- lapply(model_specs, function(spec) {
  fit_tte_model(
    data = tte_data,
    curveid = if (is.null(spec$curveid)) NULL else eval(spec$curveid, envir = tte_data),
    separate = spec$separate,
    upperl = spec$upperl
  )
})

## Access summaries: "$overall" "$treatment" "$region" and "$region_treatment"
summary(model_fits$region_treatment, robust = TRUE, units = Units)


# Loop over model list and save diagnostics to .txt files
for (model_name in names(model_fits)) {
  model <- model_fits[[model_name]]
  
  # Set up the output file path using here()
  output_txt_path <- here("Results", paste0("tte_diagnostics_", model_name, ".txt"))
  
  # Start capturing output to text file
  sink(output_txt_path)
  
  cat("=== Model Summary: ", model_name, " ===\n\n")
  print(summary(model, robust = TRUE, units = Units))
  
  # Stop capturing output
  sink()
}

#### ZIBB Modelling ####
## NB: In ZIBB modelling, petri dishes replicates (PDR) (or plates) are treated as 
## independent observations, so use average germ / dish
## Count data should be treated as numeric

## set Factors and make control the reference treatment: 
plate_summary <- plate_summary %>%
  mutate(across(c(Treatment, Pretreatment, Region, Pop, Fam), as.factor)) %>%
  mutate(Treatment = relevel(Treatment, ref = "Control"))

## double check everything has worked so far: 
str(plate_summary)

## final ZIM model: zero-inflated beta-binomial GLMM:
## proportion of germination per repliacte 
final<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Treatment * Region + Pretreatment+XR_Viability + (1 | Region / Pop),
  ziformula = ~ XR_Viability ,   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)
summary(final)

## test against other ZIBB models

#simplify models: 
# t1: drop interaction
t1<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Treatment+ Region + Pretreatment+XR_Viability + (1 | Region / Pop),
  ziformula = ~ XR_Viability ,   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)

#t2: drop  pretreatment: 
t2<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Treatment * Region +XR_Viability + (1 | Region / Pop),
  ziformula = ~ XR_Viability ,   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)

#t3: drop XR_V from conditional:
t3<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Treatment * Region + Pretreatment+ (1 | Region / Pop),
  ziformula = ~ XR_Viability ,   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)

#t4-7, complexity:
#t4: add Pretreatment* Treatment interaction: 
t4<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Pretreatment* Treatment+ Region+XR_Viability+(1 | Region / Pop),
  ziformula = ~ XR_Viability,   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)

#t5: three-way interaction
t5<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Pretreatment*Treatment*Region + XR_Viability + (1 | Region / Pop),
  ziformula = ~XR_Viability ,   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)

#t6: Pretreatments in ZI model,
t6<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Pretreatment*Treatment + Region + XR_Viability + (1 | Region / Pop),
  ziformula = ~Pretreatment + XR_Viability ,   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)

#t7 Pretreatment and Treatment in ZI
t7<-glmmTMB(
  cbind(total_germinated, Starting_number - total_germinated) ~ 
    Pretreatment*Treatment + Region + XR_Viability + (1 | Region / Pop),
  ziformula = ~Pretreatment + Treatment + XR_Viability +  (1 | Region / Pop),   # optional predictors of structural zero status
  family = betabinomial(link = "logit"),
  data = plate_summary)


BIC(final,t1,t2,t3,t4,t5,t6,t7)


## test against binomial GLMM
mod_binom <- glmmTMB(
  cbind(total_germinated,Starting_number-total_germinated)~Treatment*Region+Pretreatment +XR_Viability+(1|Region/Pop),
  family = betabinomial(link = "logit"), data = plate_summary)

# Compare AIC/BIC 
aic_comp <- AIC(final, mod_binom)
bic_comp <- BIC(final, mod_binom)


# DHARMa residual diagnostics
res_final <- simulateResiduals(final)
plot(res_final, rank = TRUE)

# Save DHARMa plots
png(filename = here("Results", "DHARMa_residuals_final_ZIM_Model.png"), width = 800, height = 600)
plot(res_final, rank = TRUE)
dev.off()

### Function to save model outputs and diagnostic tests using sink() and here(): 
save_model_output <- function(model, model_name, output_file) {
  sink(here("Results", output_file))
  
  cat(paste0("=== ", model_name, " ===\n\n"))
  print(summary(model))
  
  cat("\n--- DHARMa Residual Tests ---\n")
  res <- simulateResiduals(model)
  print(testUniformity(res))
  print(testDispersion(res))
  print(testZeroInflation(res))
  
  sink()  # Close the connection
}

# Apply the function to both models
save_model_output(final, "Final Zero-Inflated Beta-Binomial Model", "final_ZIBB_model_summary.txt")
save_model_output(mod_binom, "Binomial model without ZI", "binomial_model_summary.txt")