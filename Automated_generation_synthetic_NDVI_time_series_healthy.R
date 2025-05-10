###############################################################################
# Title:       Synthetic NDVI Time Series Generator for healthy trees (without Bark Beetle Simulation)
# Description: This script generates synthetic NDVI time series data for multiple 
#              scenarios using variable noise, trend, amplitude, seasonal phase 
#              shifts, and a simulated bark beetle infestation event. 
#              It also introduces random NA values and ensures minimum NDVI thresholds.
# Approach:    The script is a Fully Factorial Experimental Design that uses a combination of sine
#              functions to simulate seasonal variations in NDVI values, with added noise,
#              variations, linear trends, and random noise. There is not a progressive decline in
#               NDVI The script generates a CSV file for each combination of parameters

#
# Author:      Paul Arellano
# Date:        2025-04-27
#
# Inputs:      - Ranges for noise, trend, post-event decline, NA %, amplitude, etc.
#              - Time range: 2019-01-01 to 2025-12-31
#
# Outputs:     - CSV files for each combination of parameters
#              - Optionally: plots showing NDVI trends with outbreak markers
#
# Dependencies:
#              - tidyverse
#              - lubridate
#
# Usage Notes:
#              - The offset added to the seasonal signal is randomly selected 
#                between 0.01 and 0.1.
#              - Bark beetle event is simulated starting on 2022-05-15.
#              - Missing values are introduced and replaced by zeros.
###############################################################################


# Load libraries
library(tidyverse)
library(lubridate)
library(scales)
library(ggplot2)
library(mcp)
library(tibble)
library(dplyr)
library(EnvCpt)
library(segmented)


# This script generates synthetic NDVI data without a Bark Beetle Infestation outbreaks to simulate a time series for testing change-point detection models.

set.seed(42)  # Ensures reproducibility

# Time sequence: daily from Jan 2019 to Dec 2025
dates <- seq(ymd("2019-01-01"), ymd("2025-12-31"), by = "day")
n <- length(dates)

# Script Synthetic NDVI Time Series Generator without Bark Beetle Simulation ===========================================

# Define variable ranges
noise_range <- seq(from = 0.05, to = 0.15, by = 0.05)
trend_range <- seq(from = -0.05, to = 0.05, by = 0.05)
#base_range <- seq(from = 0.3, to = 0.5, by = 0.1)
#post_event_decline_range <- seq(from = -0.2, to = -0.4, by = -0.1)
random_NA_range <- seq(from = 0.05, to = 0.15, by = 0.05)
amplitude_variability <- seq(from = 0.3, to = 0.6, by = 0.1)
phase_shift_range <- c(-0.05, -0.02, 0.02, 0.05)  # Predefined phase shifts

# Time sequence: daily from Jan 2019 to Dec 2025
dates <- seq(ymd("2019-01-01"), ymd("2025-12-31"), by = "day")
n <- length(dates)

# Generate combinations of variables 
combinations <- expand.grid(
  noise = noise_range,
  trend = trend_range,
  #  base = base_range,
  #  post_event_decline = post_event_decline_range,
  random_NA = random_NA_range,
  amplitude = amplitude_variability,
  phase_shift = phase_shift_range
)

# Ensure combinations are unique
combinations <- unique(combinations)
print(paste("Final unique combinations:", nrow(combinations)))

setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Healthy_dataset")  

# Loop through each combination and generate time series
for (i in 1:nrow(combinations)) {
  # Extract current combination of variables
  noise_sd <- combinations$noise[i]
  trend_rate <- combinations$trend[i]
  #  base_value <- combinations$base[i]
  #  post_event_decline <- combinations$post_event_decline[i]
  na_percent <- combinations$random_NA[i]
  amplitude <- combinations$amplitude[i]
  phase_shift <- combinations$phase_shift[i]  # Fixed phase shift for this combination
  
  # Generate trend
  trend <- seq(from = 0, to = trend_rate, length.out = n)
  
  # Generate noise
  noise <- rnorm(n, mean = 0, sd = noise_sd)
  
  # Generate seasonal cycle with amplitude and fixed phase shift
  # Seasonal_cycle <- amplitude * sin(2 * pi * ((1:n) / 365 + phase_shift))
  # seasonal_cycle <- 0.1 + amplitude * sin(2 * pi * ((1:n) / 365 + phase_shift))
  seasonal_cycle <- amplitude * sin(2 * pi * ((1:n) / 365 + phase_shift))
  # Set a random offset between 0.02 and 0.06
  offset <- runif(1, min = 0.001, max = 0.002)
  
  # Shift seasonal cycle so its minimum becomes the random offset
  min_seasonal <- min(seasonal_cycle)
  seasonal_cycle <- seasonal_cycle + (offset - min_seasonal)
  
  # Combine base NDVI value, seasonal cycle, trend, and noise
  ndvi <- 0.3 * seasonal_cycle + trend + noise
  
  # Apply minimum threshold for NDVI values (except for NAs)
  ndvi <- pmax(ndvi, offset)
  
  # Bark beetle event: create a progressive NDVI drop
  #event_date <- ymd("2022-05-15")
  #event_index <- which(dates >= event_date)
  #decline <- rep(0, n)
  #decline[event_index] <- seq(0, post_event_decline, length.out = length(event_index))
  #ndvi <- ndvi + decline
  
  # Add missing data (NA)
  set.seed(42)  # For reproducibility
  na_count <- floor(na_percent * n)
  na_indices <- sample(1:n, na_count)
  ndvi[na_indices] <- NA
  
  # Replace NA values with 0
  ndvi[is.na(ndvi)] <- 0
  
  # Clip final NDVI values to realistic range [0, 1]
  # ndvi <- pmax(pmin(ndvi, 1), 0) # Ensure NDVI values are between 0 and 1 by clamping to 1 and 0.
  ndvi <- ifelse(ndvi >= 0 & ndvi <= 1, ndvi, NA)
  # ndvi[na_indices] <- 0  # Ensure NAs remain 0
  
  # Build dataframe
  ndvi_df <- tibble(
    date = dates,
    year = year(dates),
    day_of_year = yday(dates),
    ndvi = ndvi,
    day = seq_len(n),  # Replacing row_number() with seq_len(n)
    
    # Repeat the combination values for each row explicitly
    noise_sd = rep(noise_sd, n),        # Repeat the current noise_sd for all rows
    trend_rate = rep(trend_rate, n),    # Repeat the current trend_rate for all rows
    na_percent = rep(na_percent, n),    # Repeat the current na_percent for all rows
    amplitude = rep(amplitude, n),      # Repeat the current amplitude for all rows
    phase_shift = rep(phase_shift, n)   # Repeat the current phase_shift for all rows
  )
  
  
  # Define output file name based on variable values
  file_name <- paste0(
    "NDVI_noise", noise_sd,
    "_trend", trend_rate,
    #    "_base", base_value,
    #    "_decline", abs(post_event_decline),
    "_NA", na_percent,
    "_amplitude", amplitude,
    "_phase", phase_shift,
    ".csv"
  )
  
  #  setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Healthy_dataset")  
  
  # Save dataframe to CSV file
  write_csv(ndvi_df, file_name)
}

# Plotting time series 

# Load the generated CSV file
setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Healthy_dataset") # Replace it to your working directory 

# Prompt user to choose a CSV file from the working directory
file_name <- file.choose()  # Opens a dialog for the user to select a file
ndvi_data <- read_csv(file_name)  # Read the selected file into a dataframe
# Replace all 0 values with NA
ndvi_data[ndvi_data == 0] <- NA

# Extract the outbreak day from the dataframe
outbreak_date <- as.Date("2022-05-15")
outbreak_day <- ndvi_data %>%
  filter(date == outbreak_date) %>%
  pull(day)

x = as.Date(min(ndvi_data$date)) + 200  # Start date + offset

# Plot daily NDVI time series with outbreak annotation and file name
ggplot(ndvi_data, aes(x = date, y = ndvi)) +
  geom_line(color = "forestgreen", na.rm = TRUE) +  # Plot the NDVI time series
  scale_x_date(
    name = "Date",
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(
    name = "NDVI",
    limits = c(0,0.6)  # Set Y-scale from 0 to 6
  ) +
  labs(
    title = paste("Synthetic Daily NDVI (2019–2025) no stress events\nFile:", basename(file_name)),
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal()
