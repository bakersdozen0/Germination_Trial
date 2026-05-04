### Code to generate graphs of: 
### Figure 1) geographical distribution of sampled populations and regional assignments:
###   uses pop_coords csv, 
### Figure 5)cumulative germination over time, plotted both overall regions, 
###   and for each combination of region and treatment. 
###   Uses PCs_long.csv generated in "data transformation" section of "TTE_and_ZIBB.R" 
### Appendix 2 Scatterplots (viability metrics vs. germination)
### Bonus plots : Violin plots to display CV and variability
###             : Plot of preliminary CA trial results
### Created 25/6/2025 by JB ## 


## Library #### 
library(here)
## Data management and visualization tools:
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
## Geospatial tools: 
library(sf)
library(ncdf4)
library(raster)


#### Load required data: ####
## Load csv from parent directory if not loaded already: 
PCs_long<-read.csv(here("Proc_data","PCs_long.csv"),row.names = NULL)
pop_coors<-read.csv(here("Raw_data","Pop_coors.csv"))


##### Figure 1: map of pops / regions #### 
uk_map <- map_data("world", region = "UK")
uk_map <- uk_map %>% 
  filter(subregion != "Northern Ireland")

pop_map<-ggplot(data = uk_map) +
  geom_polygon(aes(x=long,y=lat,group=group),fill = "gray", color = "black") +  # Map of the UK with light blue fill
  geom_point(data = pop_coors, aes(x = Lon, y = Lat, color = Region), size = 3) +  # Points for populations
  geom_text_repel(data = pop_coors,aes(x = Lon, y = Lat, label = Site),vjust = -1.5, size = 3) +  
  theme_minimal() +
  #labs(title="Map of sampled populations")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill=NA,linewidth = 1),
        panel.background = element_rect(fill = "gray90", color = NA)) +
  coord_fixed(ratio = 1)  # Fix aspect ratio to ensure map is not squished

# save output in results folder: 
ggsave(
  filename = here("Results", "uk_population_map.png"),
  plot = pop_map,
  width = 8, height = 6, dpi = 300,  bg = "white" )

## Figure 4: Functions to calculate cumulative germination over time from wide-format data: ####
# Vector of regions to process
regions <- c("Lake District", "S. England", "Highlands", "S. Uplands")
## make sure start_date is saved as date, then make _weeks_end column:

PCs_long$Stratification_start_date <- lubridate::ymd(PCs_long$Stratification_start_date)
PCs_long$Check_Date <- lubridate::ymd(PCs_long$Check_Date)
#PCs_long$end_stratification <- lubridate::ymd(PCs_long$end_stratification)

PCs_long <- PCs_long %>%
  mutate(standardized_weeks_end = as.numeric(difftime(Check_Date, Stratification_start_date, units = "weeks"))) %>%
  filter(!is.na(standardized_weeks_end))  # Remove any date issues

# Function to calculate cumulative germination for a given region
calc_cum_germ <- function(region_name) {
  PCs_long %>%
    { if(region_name != "All") filter(., Region == region_name) else . } %>% 
    group_by(standardized_weeks_end, Treatment) %>%  
    summarise(Total_Germinated = sum(Germinated, na.rm = TRUE), .groups = "drop") %>% 
    arrange(standardized_weeks_end) %>% 
    group_by(Treatment) %>% 
    mutate(Cumulative_Germinated = cumsum(Total_Germinated)) %>% 
    ungroup() %>%
    mutate(Region = region_name)  # add region column for plotting
}

# Calculate for all data (no filter) + each region
all_cumulative <- calc_cum_germ("All")

regional_cumulative <- map_df(regions, calc_cum_germ)  # combined into one df

# Transition lines data frame (unchanged)
transition_lines <- data.frame(
  Transition = c("C","W->C","W->C"),
  Transition_week = c(15, 8, 23),
  Line_type = c("End of cold stratification","End of warm stratification","End of Warm, then cold stratification")
)

# Plot total cumulative germination over time (all regions combined)
totalcumgermplot <- ggplot(all_cumulative, aes(x = standardized_weeks_end, y = Cumulative_Germinated, color = Treatment)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Cumulative germination over time over all treatments and regions",
    x = "Weeks since start of stratification",
    y = "Cumulative Germinated Seeds"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("W" = "orange", "W->C" = "red", "C" = "blue", "Control" = "black")) +
  geom_vline(data = transition_lines, aes(xintercept = Transition_week, linetype = Line_type, color = Line_type), show.legend = TRUE, size = 1) +
  scale_linetype_manual(
    values = c(
      "End of cold stratification" = "dotdash", 
      "End of warm stratification" = "dashed", 
      "End of warm, then cold stratification" = "dotted"
    ),
    breaks = c(
      "End of warm stratification", 
      "End of cold stratification", 
      "End of warm, then cold stratification"
    )
  ) +
  xlim(0, 80) +
  ylim(0, 1200) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +  
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5)

print(totalcumgermplot)
ggsave(filename=here("Results","overall_germination_TS_0-80wk.png"), 
       plot=totalcumgermplot,width = 8, height = 6, dpi = 300,bg = "white")

# Plotting function for regional plots
plot_region_germ <- function(df_region) {
  region_name <- unique(df_region$Region)
  
  ggplot(df_region, aes(x = standardized_weeks_end, y = Cumulative_Germinated, color = Treatment)) +
    geom_line() +
    geom_point() +
    labs(title = region_name) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_color_manual(values = c("W" = "orange", "W->C" = "red", "C" = "blue", "Control" = "black")) +
    geom_vline(data = transition_lines, aes(xintercept = Transition_week, linetype = Line_type, color = Line_type), show.legend = TRUE, size = 1) +
    scale_linetype_manual(
      values = c(
        "End of cold stratification" = "dotdash", 
        "End of warm stratification" = "dashed", 
        "End of warm, then cold stratification" = "dotted"
      )
    ) +
    xlim(0, 80) +
    ylim(0, 400) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +  
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5)
}



# Ensure Regions are ordered alphabetically
regional_cumulative$Region <- factor(regional_cumulative$Region, levels = sort(unique(regional_cumulative$Region)))

# Create a list of plots, one per region
region_plots <- regional_cumulative %>%
  split(.$Region) %>%
  map(plot_region_germ)

regional_combined_plot <- wrap_plots(region_plots, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Cumulative germination over time for each region and treatment",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 28))
  ) & 
  labs(
    x = "Weeks since start of stratification",
    y = "Cumulative germinated seeds"
  ) & 
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10, angle = 90),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    plot.background = element_rect(fill = NA, color = NA),
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid.major = element_line(color = "gray80")
  )

print(regional_combined_plot)

ggsave(filename=here("Results","regional_germination_TS_0-80wk.png"), 
       plot=regional_combined_plot,width = 15, height = 10, units="in",dpi = 300,bg = "white")


### Appendix 2: Scatter plots of seed lot quality/ viability metrics vs germination by fam: ####
## use function to plot all required comparisons: 
create_viability_scatterplot <- function(data, x_var, y_var, x_lab, y_lab, title_text, x_lim=c(0,1), y_lim=c(0,1)) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], shape = Region, label = Fam)) +
    geom_point(aes(color = Pop), size = 2.5, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +
    geom_text_repel(aes(color = Pop), max.overlaps = 5, size = 3) + # Limit overlaps
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(x = x_lab, y = y_lab, title = title_text) +
    theme_minimal() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
}

## define plots to create with this function: 
plot_definitions <- list(
  list(file="germ_vs_fct_filled.png", x="FCT_percent", y="germ_percent", xlab="FCT Proportion Filled", ylab="Family Germination %", title="Germination vs. FCT Filled"),
  list(file="germ_vs_fct_green.png", x="FCT_green_berries", y="germ_percent", xlab="FCT Green Berries Filled", ylab="Family Germination %", title="Germination vs. FCT Green Berries"),
  list(file="germ_vs_filled_frac.png", x="Filled_fraction", y="germ_percent", xlab="Proportion Filled Seeds", ylab="Family Germination %", title="Germination vs. Filled Fraction", y_lim=c(0,0.55)),
  list(file="fct_vs_filled_frac.png", x="FCT_percent", y="Filled_fraction", xlab="FCT Proportion Filled", ylab="Filled Seed Weight Ratio", title="FCT vs. Filled Fraction", y_lim=c(0,4)),
  list(file="germ_vs_tz.png", x="TZ_Viability", y="germ_percent", xlab="TZ Viability", ylab="Family Germination %", title="Germination vs. TZ Viability", y_lim=c(0,0.55)),
  list(file="germ_vs_xr.png", x="XR_Viability", y="germ_percent", xlab="X-Ray Viability", ylab="Family Germination %", title="Germination vs. X-Ray Viability", y_lim=c(0,0.55)),
  list(file="tz_vs_xr.png", x="XR_Viability", y="TZ_Viability", xlab="X-Ray Viability", ylab="TZ Viability", title="TZ vs. X-Ray Viability", y_lim=c(0,0.8))
)

# --- Loop to create and save each plot ---
for (p in plot_definitions) {
  current_plot <- create_viability_scatterplot(
    data = full_fam_germ,
    x_var = p$x, y_var = p$y,
    x_lab = p$xlab, y_lab = p$ylab,
    title_text = p$title,
    x_lim = if(!is.null(p$x_lim)) p$x_lim else c(0,1),
    y_lim = if(!is.null(p$y_lim)) p$y_lim else c(0,1)
  )
  ggsave(
    filename = here("Results", p$file),
    plot = current_plot,
    width = 8, height = 6, dpi = 300, bg = "white"
  )
}

#### Bonus plot: Coefficient of variance visualized as violin plots: ####
## note very high variance and suspected zero-inflation (ZI is relative to model, not data) 
### Very large range of CV values by Population (13-73%): Quick visualization: 
pop_cv <- plate_summary %>%
  group_by(Pop) %>%
  summarise(cv_germ = sd(Germ_Percent)/mean(Germ_Percent)*100) %>%
  arrange(desc(cv_germ))

pop_cv_plot <- ggplot(plate_summary, aes(x = reorder(Pop, Germ_Percent, FUN = function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)), y = Germ_Percent)) +
  geom_violin(fill = "lightblue", trim = FALSE) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5, color = "darkblue") + # height = 0 keeps points on the line
  coord_flip() +
  labs(
    y = "Germination Percentage (%)", # y-axis is now Germ_Percent
    x = "Population (Ordered by CV)", # x-axis is Pop
    title = "Distribution of Germination Percentage by Population",
    subtitle = "Populations are ordered by their coefficient of variation (CV)"
  ) +
  theme_minimal()

print(pop_cv_plot)


#### Bonus plot 2; preliminary CA trial results ####
ca_data <- read.csv(here("Raw_data", "CA_exp.csv"))

# 4. Calculates the total germination and final percentage for each plate.
ca_summary <- ca_data %>%
  # Pivot all columns that start with "G." or "M."
  pivot_longer(
    cols = starts_with(c("G.", "M.")),
    names_to = "measurement_date",
    values_to = "count",
    values_drop_na = TRUE
  ) %>%
  # Separate the prefix (G/M) from the date string
  separate(measurement_date, into = c("measurement_type", "date_string"), sep = "\\.", extra = "merge") %>%
  # We only care about germination for this plot
  filter(measurement_type == "G") %>%
  # Group by each plate to sum up all germination events over time
  group_by(Region, Pop, Fam, CA_treat, Plate, Plate_ID, Starting_no) %>%
  summarise(
    total_germinated = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Calculate the final germination percentage
  mutate(
    germ_percent = total_germinated / Starting_no
  )

# --- Create the Plot ---
# This pipeline calculates the final germination percentage per plate
ca_summary <- ca_data %>%
  pivot_longer(
    cols = starts_with(c("G.", "M.")),
    names_to = "measurement_date",
    values_to = "count",
    values_drop_na = TRUE
  ) %>%
  separate(measurement_date, into = c("measurement_type", "date"), sep = "\\.", extra = "merge") %>%
  filter(measurement_type == "G") %>%
  group_by(Region, Pop, Fam, CA_treat, Plate, Plate_ID, Starting_no) %>%
  summarise(total_germinated = sum(count, na.rm = TRUE), .groups = "drop") %>%
  mutate(germ_percent = total_germinated / Starting_no)

# --- Data Augmentation for Plotting ---
# 1. Calculate the number of unique families within each population
# 2. Create a new, more informative label for the x-axis
ca_plot_data <- ca_summary %>%
  group_by(Pop) %>%
  # Add a column with the count of distinct families for each Pop
  mutate(fam_count = n_distinct(Fam)) %>%
  ungroup() %>%
  # Create a new label combining Population and family count (e.g., "BT (n=5)")
  mutate(pop_label = paste0(Pop, " (n=", fam_count, ")"))


# --- Create the Single-Panel Faceted Boxplot ---

ca_single_panel_plot <- ggplot(ca_plot_data, aes(x = reorder(pop_label, germ_percent, FUN = median), y = germ_percent)) +
  # Jittered points are plotted first to sit behind the boxplots
  geom_jitter(width = 0.25, alpha = 0.4, height = 0) +
  # Boxplots summarize the distribution (the box will be squashed at zero)
  geom_boxplot(aes(fill = Pop), alpha = 0.7, outlier.shape = NA) + # outlier.shape=NA avoids plotting outliers twice
  # Create vertical panels for each CA treatment level
  facet_grid(CA_treat ~ ., labeller = labeller(CA_treat = ~ paste("CA Conc:", .x, "%"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_flip() + # Flip coordinates to make population labels horizontal and easy to read
  labs(
    title = "Germination Distribution by Citric Acid Concentration",
    subtitle = "Populations are ordered by median germination; labels show family count (n)",
    x = "Population",
    y = "Final Germination Rate"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text.y = element_text(face = "bold", angle = 0), # Horizontal facet labels
    strip.background = element_rect(fill = "gray90"),
    axis.text.y = element_text(face = "bold") # Bold population labels
  )


## option 2: 
ca_family_trends_plot <- ggplot(ca_summary, aes(x = factor(CA_treat), y = germ_percent, group = Fam)) +
  # Add a faint line for each individual family to show variation
  stat_summary(fun = "mean", geom = "line", aes(color = Pop), alpha = 0.4) +
  # Add a thick, black line for the overall average trend across all families
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color = "black", linewidth = 1.5) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Germination Trend by Citric Acid Concentration",
    subtitle = "Each faint line represents a seed family; the black line is the overall average",
    x = "Citric Acid Treatment Concentration (%)",
    y = "Mean Germination Rate"
  ) +
  theme_bw()
ca_family_trends_plot