# ============================================================
# R Script for NDVI Change-Point Detection Using MCP Model
# ============================================================
# Author: Paul Arellano, Ph.D.
# Date: 2025-04-19
# Description:
# This script generates daily synthetic NDVI data sets for both 1) bark beetle infestation and 
# Drought stress for the following time frams: 2019-01-01 to 2025-12-31.
# The script applies the following change-point models to the time series, extracts change-point estimates, and visualizes results.
# 1.  CUSUM
# 2.  MCP
# 3.  EnvCPT
# 4.  Segmented
# 5.  STRCHANGE-breakpoints
# 6.  STRCHANGE:: fstats
# 7.  STRUCCHANGE
# 8.  bcp
# 9.  ecp
# 10. Change Point
# 11. cpt
# 12. TSMCP
# 13. cpm

# Key Steps:
# 1. Load necessary libraries (dplyr, tidyr, future, mcp, ggplot2, etc.).
# 2. Read NDVI CSV files from the specified directory.
# 3. Define the MCP model structure (initial intercept, linear trend).
# 4. Preprocess NDVI values (handle zero values, interpolate missing data).
# 5. Set up parallel processing for optimized computations.
# 6. Process each dataset:
#     - Read CSV using 'fread' for speed.
#     - Validate required columns ('ndvi', 'day').
#     - Apply NDVI preprocessing.
#     - Fit the MCP model and extract change-point estimates.
# 7. Aggregate and save results in a CSV file.
# 8. Visualize results using histograms and distribution plots.
# 9. Apply filtering to remove extreme outliers in detected change-point values.
# 10. Compute statistical summaries (change-point percentages within a predefined range).
#
# Output:
# - A summary CSV file containing change-point estimates.
# - Histograms and cumulative distribution plots of detected change-points.
#
# Notes:
# - This script is designed for synthetic NDVI datasets simulating ecological stress factors.
# - Requires the 'mcp' package for Bayesian change-point modeling.
# - Adjust file paths and model parameters as needed.
#
# ============================================================

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


# This script generates synthetic NDVI data with a Bark Beetle Infestation outbreak to simulate a time series for testing change-point detection models.

# Define variable ranges
noise_range <- seq(from = 0.05, to = 0.15, by = 0.05)
trend_range <- seq(from = -0.1, to = -0.3, by = -0.1)
post_event_decline_range <- seq(from = -0.2, to = -0.4, by = -0.1)
random_NA_range <- seq(from = 0.00, to = 0.10, by = 0.05)
amplitude_variability <- seq(from = 0.3, to = 0.6, by = 0.1)
phase_shift_range <- c(-0.05, -0.02, 0.02, 0.05)  # Predefined phase shifts

# Time sequence: daily from Jan 2019 to Dec 2025
dates <- seq(ymd("2019-01-01"), ymd("2025-12-31"), by = "day")
n <- length(dates)

# Generate combinations of variables 
combinations <- expand.grid(
  noise = noise_range,
  trend = trend_range,
  post_event_decline = post_event_decline_range,
  random_NA = random_NA_range,
  amplitude = amplitude_variability,
  phase_shift = phase_shift_range
)

# Ensure combinations are unique
combinations <- unique(combinations)
print(paste("Final unique combinations:", nrow(combinations)))

setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25")  

# Loop through each combination and generate time series
for (i in 1:nrow(combinations)) {
  # Extract current combination of variables
  noise_sd <- combinations$noise[i]
  trend_rate <- combinations$trend[i]
  post_event_decline <- combinations$post_event_decline[i]
  na_percent <- combinations$random_NA[i]
  amplitude <- combinations$amplitude[i]
  phase_shift <- combinations$phase_shift[i]  # Fixed phase shift for this combination
  
  # Generate trend
  trend <- seq(from = 0, to = trend_rate, length.out = n)
  
  # Generate noise
  noise <- rnorm(n, mean = 0, sd = noise_sd)
  
  # Seasonal cycle with amplitude and fixed phase shift
  seasonal_cycle <- amplitude * sin(2 * pi * ((1:n) / 365 + phase_shift))
  
  # Define seasonal reduction for spring & summer (March–August)
  seasonal_reduction <- ifelse(month(dates) %in% 3:8, runif(n, min = 0.02, max = 0.06), 0)
  
  # Set a random offset between 0.1 and 0.2 for fall/winter months
  offset <- ifelse(month(dates) %in% c(9:12, 1:2), runif(n, min = 0.1, max = 0.2), runif(n, min = 0.02, max = 0.06))
  
  # Adjust seasonal cycle so fall/winter values stay above 0.1
  min_seasonal <- min(seasonal_cycle)
  seasonal_cycle <- seasonal_cycle + (offset - min_seasonal)
  
  # Apply seasonal reduction
  adjusted_amplitude <- amplitude * (1 - seasonal_reduction)
  
  # Calculate final NDVI values
  ndvi <- 0.4 * seasonal_cycle + trend + noise
  
  # Ensure fall/winter NDVI stays in range [0.2, 0.3]
  ndvi <- ifelse(month(dates) %in% c(9:12, 1:2), pmax(ndvi, offset), ndvi)
  
  # Bark beetle event: create a progressive NDVI drop
  event_date <- ymd("2022-05-15")
  event_index <- which(dates >= event_date)
  decline <- rep(0, n)
  decline[event_index] <- seq(0, post_event_decline, length.out = length(event_index))
  ndvi <- ndvi + decline
  
  # Add missing data (NA)
  set.seed(42)  # For reproducibility
  na_count <- floor(na_percent * n)
  na_indices <- sample(1:n, na_count)
  ndvi[na_indices] <- NA
  
  # Replace NA values with 0
  ndvi[is.na(ndvi)] <- 0
  # Replace values less than 0 with NA
  ndvi[ndvi < 0] <- 0
  
  # Build dataframe
  ndvi_df <- tibble(
    date = dates,
    year = year(dates),
    day_of_year = yday(dates),
    ndvi = ndvi,
    day = seq_len(n)  # Replacing row_number() with seq_len(n)
  )
  
  # Define output file name based on variable values
  file_name <- paste0(
    "NDVI_noise", noise_sd,
    "_trend", trend_rate,
    "_decline", abs(post_event_decline),
    "_NA", na_percent,
    "_amplitude", amplitude,
    "_phase", phase_shift,
    ".csv"
  )
  
  # Save file
  write.csv(ndvi_df, file = file_name, row.names = FALSE)
}



# Plotting time series 

setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25")

# Prompt user to choose a CSV file from the working directory
file_name <- file.choose()  # Opens a dialog for the user to select a file
ndvi_data <- read_csv(file_name)  # Read the selected file into a dataframe

# Extract the outbreak day from the dataframe
outbreak_date <- as.Date("2022-05-15")
outbreak_day <- ndvi_data %>%
  filter(date == outbreak_date) %>%
  pull(day)

x = as.Date(min(ndvi_data$date)) + 200  # Start date + offset

# Plot daily NDVI time series with outbreak annotation and file name
ggplot(ndvi_data, aes(x = date, y = ndvi)) +
  geom_line(color = "forestgreen", na.rm = TRUE) +  # Plot the NDVI time series
  geom_vline(xintercept = outbreak_date, 
             color = "red", linetype = "dashed", linewidth = 1) +  # Outbreak line
  annotate("text", x = outbreak_date + 10, y = 0.9,  # Print outbreak day
           label = paste("Outbreak day:", outbreak_day), 
           color = "red", size = 4, hjust = 0) +
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
    title = paste("Synthetic Daily NDVI (2019–2025) with Bark Beetle Infestation\nFile:", basename(file_name)),
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal()


# Plot a subset of the daaset

# Filter data for the selected time range
ndvi_subset <- ndvi_data %>%
  filter(date >= ymd("2020-01-01") & date <= ymd("2021-09-30"))

# Create the plot
ggplot(ndvi_subset, aes(x = date, y = ndvi)) +
  geom_line(color = "darkgreen", size = 1) +  # Line plot
#  geom_point(color = "black", size = 2, alpha = 0.7) +  # Add points for visibility
  labs(
    title = "NDVI Time Series (Jan 2020 - Sep 2021)",
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
# ==============================================================================
#
#
#
#
# ==============================================================================
# Synthetic NDVI data with a drought event to simulate a time series for testing change-point detection models.

# Time sequence: weekly from Jan 2019 to Dec 2025
dates <- seq(ymd("2019-01-01"), ymd("2025-12-31"), by = "day")
n <- length(dates)

# Seasonal NDVI pattern: 1 cycle per year
seasonal_cycle <- sin(2 * pi * (1:n) / 365)

# Long-term slight trend: Gradual vegetation loss
trend <- seq(0, -0.4, length.out = n)

# Random noise
noise <- rnorm(n, mean = 0, sd = 0.05)

# Drought event: 2-year decline
drought_start_date <- ymd("2021-06-01")
drought_end_date <- ymd("2023-06-01")
drought_index <- which(dates >= drought_start_date & dates <= drought_end_date)

# Recovery period: gradual recovery for 3 years after the drought
recovery_start_date <- ymd("2023-06-01")
recovery_end_date <- ymd("2026-06-01")  # 3 years of recovery
recovery_index <- which(dates > recovery_start_date & dates <= recovery_end_date)

# Gradual NDVI drop during drought (allowing low values but avoiding 0)
decline <- rep(0, n)
decline[drought_index] <- seq(0.1, 0, length.out = length(drought_index))  # Decline starts at 0.1 and remains above 0

# Progressive seasonality decline during drought
seasonal_amplitude <- rep(1, n)
seasonal_amplitude[drought_index] <- seq(1, 0.1, length.out = length(drought_index))  # Minimum amplitude is 0.1

# Gradual recovery in seasonal amplitude post-drought
seasonal_amplitude[recovery_index] <- seq(0.1, 1, length.out = length(recovery_index))  # Recovery from 0.1 to full amplitude

# Recovery in overall NDVI values post-drought
recovery <- rep(0, n)
recovery[recovery_index] <- seq(0, 0.5, length.out = length(recovery_index))  # Gradual recovery from low to normal

# Apply the adjustments
adjusted_seasonal_cycle <- seasonal_amplitude * seasonal_cycle
ndvi <- 0.5 + 0.3 * adjusted_seasonal_cycle + trend + decline + recovery + noise

# Clip NDVI values to realistic range [0.1, 1] (ensuring a minimum of 0.1)
ndvi <- pmax(pmin(ndvi, 1), 0.1)

# Build dataframe
ndvi_df_drought <- tibble(
  date = dates,
  year = year(dates),
  week = isoweek(dates),
  ndvi = ndvi
)

#Creates a sequential time variable
ndvi_df_drought <- ndvi_df_drought %>%
  mutate(day = row_number())

# Adjust recovery start date to begin after the drought end date
recovery_start_date <- drought_end_date + days(1)

# Extract the day values for drought start and recovery start dates
drought_start_day <- ndvi_df_drought %>%
  filter(date == drought_start_date) %>%
  pull(day)

recovery_start_day <- ndvi_df_drought %>%
  filter(date == recovery_start_date) %>%
  pull(day)

# Plot the NDVI time series with annotations
ggplot(ndvi_df_drought, aes(x = date, y = ndvi)) +
  geom_line(color = "forestgreen") +
  geom_vline(xintercept = drought_start_date, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = recovery_start_date, color = "green", linetype = "dashed", linewidth = 1) +  # Green line shifted
  annotate("text", x = drought_start_date + 70, y = 0.9,
           label = paste("Start of Drought\nDay:", drought_start_day), color = "red", size = 4, hjust = 0) +
  annotate("text", x = recovery_start_date + 70, y = 0.9,
           label = paste("Recovery Begins\nDay:", recovery_start_day), color = "green", size = 4, hjust = 0) +
  scale_x_date(
    name = "Date",
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  labs(
    title = "Synthetic NDVI with Two-Year Decline and Recovery",
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal()

# Save to CSV
write.csv(ndvi_df_drought, "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/MCP_model/Synthetic_data_NDVI_drought_event_2.csv", row.names = FALSE)
# ==================================================================================================
#
#
#
# ====================== CUSUM MODEL Beetle ========================================================

# Ensure there are no NA values in the ndvi column for calculations
ndvi_clean <- ndvi_df_beetle$ndvi
ndvi_clean[is.na(ndvi_clean)] <- mean(ndvi_clean, na.rm = TRUE)  # Replace NA with column mean

# Calculate CUSUM values
cusum_vals <- cumsum(ndvi_clean - mean(ndvi_clean, na.rm = TRUE))

# Calculate positive changes
cusum_positive <- cusum_vals[cusum_vals > 0 & !is.na(cusum_vals)]  # Exclude NA values
total_positive_change <- sum(cusum_positive, na.rm = TRUE)  # Sum of positive changes

# Calculate negative changes
cusum_negative <- cusum_vals[cusum_vals < 0 & !is.na(cusum_vals)]  # Exclude NA values
total_negative_change <- sum(cusum_negative, na.rm = TRUE)  # Sum of negative changes

# Identify the change point (index of the largest absolute change in CUSUM)
change_point <- which.max(abs(cusum_vals))

# Plot the CUSUM values
plot(cusum_vals, type = "l", col = "blue", main = "CUSUM Values", xlab = "Index", ylab = "CUSUM")
# Add a vertical dashed line at the change point
abline(v = change_point, col = "red", lty = 2)
# Highlight the change point with a red dot
points(change_point, cusum_vals[change_point], col = "red", pch = 19)
# Add a label next to the red point
text(x = change_point + 10,  # Adjust position slightly to the right of the point
     y = cusum_vals[change_point],  # Align with the y-value of the point
     labels = paste0("Change Point: ", change_point), 
     col = "red", pos = 4, cex = 0.8)  # pos = 4 places the text to the right

# Print the results
print(paste("Total Positive Change:", total_positive_change))
print(paste("Total Negative Change:", total_negative_change))
print(paste("Change Point Index:", change_point))

# Plotting
ggplot() +
  # Scatter points for NDVI
  geom_point(aes(x = ndvi_df_beetle$weeks, y = ndvi_df_beetle$ndvi), color = "black") +
  # Smooth blue line for NDVI trend
  geom_smooth(aes(x = ndvi_df_beetle$weeks, y = ndvi_df_beetle$ndvi), method = "loess", color = "blue", linewidth = 1) +
  # Green area for positive CUSUM values
  geom_area(data = tibble(x = ndvi_df_beetle$weeks, y = cusum_vals),
            aes(x = x, y = ifelse(cusum_vals > 0, cusum_vals, 0)), fill = "green", alpha = 0.5) +
  # Red area for negative CUSUM values
  geom_area(data = tibble(x = ndvi_df_beetle$weeks, y = cusum_vals),
            aes(x = x, y = ifelse(cusum_vals < 0, cusum_vals, 0)), fill = "red", alpha = 0.5) +
  # Vertical dashed line for change point
  geom_vline(xintercept = ndvi_df_beetle$weeks[change_point], color = "black", linetype = "dashed", size = 1.2) +
  # Labels and axis formatting
  labs(
    title = "Synthetic Dataset - CUSUM Positive and Negative Changes",
    x = "Weeks",
    y = "CUSUM Value"
  ) +
  scale_x_continuous(
    breaks = seq(min(ndvi_df_beetle$weeks), max(ndvi_df_beetle$weeks), by = 10)
  ) +
  theme_minimal()


# ==============================================================================

# ====================== CUSUM MODEL Drought  ========================================================

# Ensure there are no NA values in the ndvi column for calculations
ndvi_clean <- ndvi_df_drought$ndvi
ndvi_clean[is.na(ndvi_clean)] <- mean(ndvi_clean, na.rm = TRUE)  # Replace NA with column mean

# Calculate CUSUM values
cusum_vals <- cumsum(ndvi_clean - mean(ndvi_clean, na.rm = TRUE))

# Calculate positive changes
cusum_positive <- cusum_vals[cusum_vals > 0 & !is.na(cusum_vals)]  # Exclude NA values
total_positive_change <- sum(cusum_positive, na.rm = TRUE)  # Sum of positive changes

# Calculate negative changes
cusum_negative <- cusum_vals[cusum_vals < 0 & !is.na(cusum_vals)]  # Exclude NA values
total_negative_change <- sum(cusum_negative, na.rm = TRUE)  # Sum of negative changes

# Identify the change point (index of the largest absolute change in CUSUM)
change_point <- which.max(abs(cusum_vals))

# Plot the CUSUM values
plot(cusum_vals, type = "l", col = "blue", main = "CUSUM Values", xlab = "Index", ylab = "CUSUM")
# Add a vertical dashed line at the change point
abline(v = change_point, col = "red", lty = 2)
# Highlight the change point with a red dot
points(change_point, cusum_vals[change_point], col = "red", pch = 19)
# Add a label next to the red point
text(x = change_point + 10,  # Adjust position slightly to the right of the point
     y = cusum_vals[change_point],  # Align with the y-value of the point
     labels = paste0("Change Point: ", change_point), 
     col = "red", pos = 4, cex = 0.8)  # pos = 4 places the text to the right

# Print the results
print(paste("Total Positive Change:", total_positive_change))
print(paste("Total Negative Change:", total_negative_change))
print(paste("Change Point Index:", change_point))

# Plotting
ggplot() +
  # Scatter points for NDVI
  geom_point(aes(x = ndvi_df_drought$weeks, y = ndvi_df_drought$ndvi), color = "black") +
  # Smooth blue line for NDVI trend
  geom_smooth(aes(x = ndvi_df_drought$weeks, y = ndvi_df_drought$ndvi), method = "loess", color = "blue", linewidth = 1) +
  # Green area for positive CUSUM values
  geom_area(data = tibble(x = ndvi_df_drought$weeks, y = cusum_vals),
            aes(x = x, y = ifelse(cusum_vals > 0, cusum_vals, 0)), fill = "green", alpha = 0.5) +
  # Red area for negative CUSUM values
  geom_area(data = tibble(x = ndvi_df_drought$weeks, y = cusum_vals),
            aes(x = x, y = ifelse(cusum_vals < 0, cusum_vals, 0)), fill = "red", alpha = 0.5) +
  # Vertical dashed line for change point
  geom_vline(xintercept = ndvi_df_drought$weeks[change_point], color = "black", linetype = "dashed", size = 1.2) +
  # Labels and axis formatting
  labs(
    title = "Synthetic Dataset - CUSUM Positive and Negative Changes",
    x = "Weeks",
    y = "CUSUM Value"
  ) +
  scale_x_continuous(
    breaks = seq(min(ndvi_df_drought$weeks), max(ndvi_df_drought$weeks), by = 10)
  ) +
  theme_minimal()

# ============================  MCP ==================================================

# Define the MCP model. GAUSSIAN. Change point model with three segments 
model <- list(
  ndvi ~ 1,  # Initial intercept 
  #~ day,    # Linear trend before the change point.
  ~ day     # Linear trend before the change point. 
)

# Fit the model to the ndvi_df dataset
fit = mcp(model, data = ndvi_df_drought)
#fit = mcp(model, data = deseasonalized_short_trend)
# Plot the fitted model
plot(fit, q_fit = TRUE) # Plot the fitted model with change points, including the 95% Highest Density Interval (HDI) which is the natural Bayesian credible interval (q_fit=TRUE)
# Set global options to avoid scientific notation
options(scipen = 999)  # Increase the penalty for scientific notation
# Summary of the fit
summary(fit)
# Plot the fitted model with change points
plot(fit, q_fit = TRUE) # Plot the fitted model with change points, including the 95% Highest Density Interval (HDI) which is the natural Bayesian credible interval (q_fit=TRUE)

# Plot the change points
plot(fit) + plot_pars(fit, pars = c("cp_1"), type = "dens_overlay")


# MCP model for the drought event ====================================

# Define the MCP model. GAUSSIAN. Change point model with three segments 
model <- list(
  ndvi ~ 1,  # Initial intercept 
  #~ day,    # Linear trend before the change point.
  ~ day     # Linear trend before the change point. 
)

# Fit the model to the ndvi_df dataset
fit = mcp(model, data = ndvi_df_beetle)
#fit = mcp(model, data = deseasonalized_short_trend)
# Plot the fitted model
plot(fit, q_fit = TRUE) # Plot the fitted model with change points, including the 95% Highest Density Interval (HDI) which is the natural Bayesian credible interval (q_fit=TRUE)
# Set global options to avoid scientific notation
options(scipen = 999)  # Increase the penalty for scientific notation
# Summary of the fit
summary(fit)
# Plot the fitted model with change points
plot(fit, q_fit = TRUE) # Plot the fitted model with change points, including the 95% Highest Density Interval (HDI) which is the natural Bayesian credible interval (q_fit=TRUE)

# Plot the change points
plot(fit) + plot_pars(fit, pars = c("cp_1"), type = "dens_overlay")

# =========================== MCP TRIGONOMETRIC Need more develop ===================================================
# Define the MCP model with trigonometric terms
model <- list(
  ndvi ~ 1,                                    # Initial intercept
  ~ day,                                       # Linear trend for segment 2
  ~ sin(2 * 3.14159 * day / 365) +             # Seasonal effect in segment 3
    cos(2 * 3.14159 * day / 365)
)

# Fit the model to the ndvi_df dataset
fit = mcp(model, data = ndvi_df_drought)
#fit = mcp(model, data = deseasonalized_short_trend)
# Plot the fitted model
plot(fit, q_fit = TRUE) # Plot the fitted model with change points, including the 95% Highest Density Interval (HDI) which is the natural Bayesian credible interval (q_fit=TRUE)
# Set global options to avoid scientific notation
options(scipen = 999)  # Increase the penalty for scientific notation
# Summary of the fit
summary(fit)
# Plot the fitted model with change points
plot(fit, q_fit = TRUE) # Plot the fitted model with change points, including the 95% Highest Density Interval (HDI) which is the natural Bayesian credible interval (q_fit=TRUE)

# Plot the change points
plot(fit) + plot_pars(fit, pars = c("cp_1"), type = "dens_overlay")


# ==================================== EnvCpt BEETLE ==================================
ndvi_vector <- ndvi_df_beetle$ndvi  # Extract the NDVI column as a numeric vector
#ndvi_vector <- deseasonalized_short_trend$trend  # Extract the NDVI column as a numeric vector
ndvi_vector_clean <- na.omit(ndvi_vector)

fit_envcpt = envcpt(data = ndvi_vector_clean)  # Fit all models at once
fit_envcpt$summary  # Show log-likelihoods
plot(fit_envcpt)

fit_envcpt$meancpt@cpts  # Show the mean change point
# Show the mean + AR2 change points
mean_ar2_cpts <- fit_envcpt$meanar2cpt@cpts
print(mean_ar2_cpts)  # Display the mean + AR2 change points
trend_ar2_cpts <- fit_envcpt$trendar2cpt@cpts
print(trend_ar2_cpts)  # Display the trend + AR2 change points
trend_ar1_cpts <- fit_envcpt$trendar1cpt@cpts
print(trend_ar1_cpts)  # Display the trend + AR2 change points
meancpt_AR2 <- fit_envcpt$meanar2cpt
print(meancpt_AR2)  # Display the mean + AR2 change points
trendcpt <- fit_envcpt$trendcpt
print(trendcpt)  # Display the mean + AR2 change points

# MEAN DIFFERENCE ==============================
mean_differences <- numeric(length(fit_envcpt$meancpt@cpts) - 1)

# Loop through change points and calculate mean differences
for (i in seq_along(fit_envcpt$meancpt@cpts)[-1]) {
  segment_1 <- ndvi_vector_clean[fit_envcpt$meancpt@cpts[i - 1]:fit_envcpt$meancpt@cpts[i]]
  segment_2 <- ndvi_vector_clean[fit_envcpt$meancpt@cpts[i]:ifelse(i == length(fit_envcpt$meancpt@cpts), length(ndvi_vector_clean), fit_envcpt$meancpt@cpts[i + 1])]
  mean_differences[i - 1] <- abs(mean(segment_1) - mean(segment_2))
}

# Identify the indices of the two largest differences
top_two_indices <- order(mean_differences, decreasing = TRUE)[1:2]

# Extract the two most significant change points
most_significant_change_points <- fit_envcpt$meancpt@cpts[top_two_indices]
print(most_significant_change_points)

# TREND AR1 ==============================

# Initialize a vector for slope differences
trend_differences <- numeric(length(trend_ar1_cpts) - 1)

# Loop through change points and calculate slope differences
for (i in seq_along(trend_ar1_cpts)[-1]) {
  segment_1 <- ndvi_vector_clean[trend_ar1_cpts[i - 1]:trend_ar1_cpts[i]]
  segment_2 <- ndvi_vector_clean[trend_ar1_cpts[i]:ifelse(i == length(trend_ar1_cpts), length(ndvi_vector_clean), trend_ar1_cpts[i + 1])]
  
  # Fit linear trends for each segment
  trend_1 <- lm(segment_1 ~ seq_along(segment_1))$coefficients[2]  # Slope of segment 1
  trend_2 <- lm(segment_2 ~ seq_along(segment_2))$coefficients[2]  # Slope of segment 2
  
  # Calculate the absolute difference in slopes
  trend_differences[i - 1] <- abs(trend_1 - trend_2)
}

# Identify the indices of the two largest differences
top_two_indices <- order(trend_differences, decreasing = TRUE)[1:2]

# Extract the two most representative change points
most_representative_changes <- trend_ar1_cpts[top_two_indices]
print(most_representative_changes)

# TREND CPT ==============================

# Initialize a vector for slope differences
trend_differences <- numeric(length(trendcpt) - 1)

# Loop through change points and calculate slope differences
for (i in seq_along(trend_cpts)[-1]) {
  segment_1 <- ndvi_vector_clean[trend_cpts[i - 1]:trend_cpts[i]]
  segment_2 <- ndvi_vector_clean[trend_cpts[i]:ifelse(i == length(trend_cpts), length(ndvi_vector_clean), trend_cpts[i + 1])]
  
  # Fit linear trends for each segment
  trend_1 <- lm(segment_1 ~ seq_along(segment_1))$coefficients[2]  # Slope of segment 1
  trend_2 <- lm(segment_2 ~ seq_along(segment_2))$coefficients[2]  # Slope of segment 2
  
  # Calculate the absolute difference in slopes
  trend_differences[i - 1] <- abs(trend_1 - trend_2)
}

# Identify the indices of the two largest differences
top_two_indices <- order(trend_differences, decreasing = TRUE)[1:2]

# Extract the two most representative change points
most_representative_changes <- trend_ar1_cpts[top_two_indices]
print(most_representative_changes)






# ======================== EnvCpt Drought ==================================

ndvi_vector_drought <- ndvi_df_drought$ndvi  # Extract the NDVI column as a numeric vector
#ndvi_vector <- deseasonalized_short_trend$trend  # Extract the NDVI column as a numeric vector
ndvi_vector_clean <- na.omit(ndvi_vector_drought)

fit_envcpt = envcpt(data = ndvi_vector_clean)  # Fit all models at once
fit_envcpt$summary  # Show log-likelihoods
plot(fit_envcpt)

fit_envcpt$meancpt@cpts  # Show the mean change point
# Show the mean + AR2 change points
mean_ar2_cpts <- fit_envcpt$meanar2cpt@cpts
print(mean_ar2_cpts)  # Display the mean + AR2 change points
trend_ar2_cpts <- fit_envcpt$trendar2cpt@cpts
print(trend_ar2_cpts)  # Display the trend + AR2 change points
trend_ar1_cpts <- fit_envcpt$trendar1cpt@cpts
print(trend_ar1_cpts)  # Display the trend + AR2 change points
meancpt_AR2 <- fit_envcpt$meanar2cpt
print(meancpt_AR2)  # Display the mean + AR2 change points


# MEAN DIFFERENCE ==============================
mean_differences <- numeric(length(fit_envcpt$meancpt@cpts) - 1)

# Loop through change points and calculate mean differences
for (i in seq_along(fit_envcpt$meancpt@cpts)[-1]) {
  segment_1 <- ndvi_vector_clean[fit_envcpt$meancpt@cpts[i - 1]:fit_envcpt$meancpt@cpts[i]]
  segment_2 <- ndvi_vector_clean[fit_envcpt$meancpt@cpts[i]:ifelse(i == length(fit_envcpt$meancpt@cpts), length(ndvi_vector_clean), fit_envcpt$meancpt@cpts[i + 1])]
  mean_differences[i - 1] <- abs(mean(segment_1) - mean(segment_2))
}

# Identify the indices of the two largest differences
top_two_indices <- order(mean_differences, decreasing = TRUE)[1:2]

# Extract the two most significant change points
most_significant_change_points <- fit_envcpt$meancpt@cpts[top_two_indices]
print(most_significant_change_points)

# TREND AR1 ==============================

# Initialize a vector for slope differences
trend_differences <- numeric(length(trend_ar1_cpts) - 1)

# Loop through change points and calculate slope differences
for (i in seq_along(trend_ar1_cpts)[-1]) {
  segment_1 <- ndvi_vector_clean[trend_ar1_cpts[i - 1]:trend_ar1_cpts[i]]
  segment_2 <- ndvi_vector_clean[trend_ar1_cpts[i]:ifelse(i == length(trend_ar1_cpts), length(ndvi_vector_clean), trend_ar1_cpts[i + 1])]
  
  # Fit linear trends for each segment
  trend_1 <- lm(segment_1 ~ seq_along(segment_1))$coefficients[2]  # Slope of segment 1
  trend_2 <- lm(segment_2 ~ seq_along(segment_2))$coefficients[2]  # Slope of segment 2
  
  # Calculate the absolute difference in slopes
  trend_differences[i - 1] <- abs(trend_1 - trend_2)
}

# Identify the indices of the two largest differences
top_two_indices <- order(trend_differences, decreasing = TRUE)[1:2]

# Extract the two most representative change points
most_representative_changes <- trend_ar1_cpts[top_two_indices]
print(most_representative_changes)

# TREND AR2 ==============================

# Initialize a vector for slope differences
trend_differences <- numeric(length(trend_ar2_cpts) - 1)

# Loop through change points and calculate slope differences
for (i in seq_along(trend_ar1_cpts)[-1]) {
  segment_1 <- ndvi_vector_clean[trend_ar1_cpts[i - 1]:trend_ar1_cpts[i]]
  segment_2 <- ndvi_vector_clean[trend_ar1_cpts[i]:ifelse(i == length(trend_ar1_cpts), length(ndvi_vector_clean), trend_ar1_cpts[i + 1])]
  
  # Fit linear trends for each segment
  trend_1 <- lm(segment_1 ~ seq_along(segment_1))$coefficients[2]  # Slope of segment 1
  trend_2 <- lm(segment_2 ~ seq_along(segment_2))$coefficients[2]  # Slope of segment 2
  
  # Calculate the absolute difference in slopes
  trend_differences[i - 1] <- abs(trend_1 - trend_2)
}

# Identify the indices of the two largest differences
top_two_indices <- order(trend_differences, decreasing = TRUE)[1:2]

# Extract the two most representative change points
most_representative_changes <- trend_ar1_cpts[top_two_indices]
print(most_representative_changes)





















# =================================================================================

# =================================SEGMENTED =========================================
# Create a new dataframe with the dependent and independent variables
data_for_lm <- data.frame(
  y = ndvi_df_beetle$ndvi,          # Dependent variable
  x = ndvi_df_beetle$day           # Independent variable (e.g., time column)
)
fit_lm <- lm(y ~ 1 + x, data = data_for_lm)  # Intercept-only model with trend (time)
fit_segmented = segmented(fit_lm, seg.Z = ~x, npsi = 1)  # Two change points along x
summary(fit_segmented)
plot(fit_segmented)
#axis(1, at = seq(0, 500, by = 5))  # Add x-axis ticks at intervals of 5
points(data_for_lm)
lines.segmented(fit_segmented)
points.segmented(fit_segmented)


# Message for confirmation
print("Change-point detection completed. Results saved as 'change_points_results.csv'.")

# ======================================== SEGMENTED drought=========================================

# Create a new dataframe with the dependent and independent variables

data_for_lm <- data.frame(
  y = ndvi_df_drought$ndvi,          # Dependent variable
  x = ndvi_df_drought$day           # Independent variable (e.g., time column)
)
fit_lm <- lm(y ~ 1 + x, data = data_for_lm)  # Intercept-only model with trend (time)

fit_segmented = segmented(fit_lm, seg.Z = ~x, npsi = 2)  # Two change points along x
summary(fit_segmented)
plot(fit_segmented)
points(data_for_lm)
lines.segmented(fit_segmented)
points.segmented(fit_segmented)

# Message for confirmation
print("Change-point detection completed. Results saved as 'change_points_results.csv'.")

# ============================================== OTHER MODELS ==========================================
ndvi_df_beetle <- ndvi_df

ndvi_df_drought <- ndvi_df

# ========================================= STRCHANGE-BREAKPOINTS MODEL =========================================

# Apply the strucchange::breakpoints function
library(strucchange)
#breakpoints_model <- breakpoints(ndvi_df_beetle$ndvi ~ 1, data = ndvi_df_beetle)
#breakpoints_model <- breakpoints(ndvi_df_drought$ndvi ~ 1, h = 100)# Set minimum segment size (h)
# Fit model with exactly 1 breakpoint
#breakpoints_model <- breakpoints(ndvi_df_beetle$ndvi ~ 1, breaks = 2)
breakpoints_model <- breakpoints(ndvi_df_beetle$ndvi ~ ndvi_df_beetle$day, breaks = 3)

# Plot the breakpoints
plot(breakpoints_model)
print(breakpoints_model)  # Print the breakpoints)

break_indices <- breakpoints_model$breakpoints  # Breakpoint positions
ndvi_at_breaks <- ndvi_df_beetle$ndvi[break_indices]  # NDVI values at breakpoints
# Corresponding day values
days_at_breaks <- ndvi_df_beetle$day[break_indices]  # Days matching the breakpoints
print(ndvi_at_breaks)
print(days_at_breaks)  # Print corresponding day values

# ================================= STRUCCHANGE::fSTATS MODEL =========================================

# This model detect only one breakpoint then it is not suitable for the drought event

library(zoo)
ndvi_df_drought$ndvi_smoothed <- rollmean(ndvi_df_drought$ndvi, k = 3, fill = NA)
# Apply the Fstats function to detect structural changes
fstats_model <- Fstats(ndvi_smoothed ~ day, data = ndvi_df_drought)
# Plot the F-statistics
plot(fstats_model, main = "F-Statistics for Structural Changes - Drought Event")
# Extract the breakpoints (2 main structural changes)
breakpoints_fstats <- breakpoints(fstats_model, breaks = 2, h = 10)  # Adjust 'h' for minimum segment size
print(breakpoints_fstats)
# Extract breakpoint indices
break_indices <- breakpoints_fstats$breakpoints
print(break_indices)
# If more than one breakpoint detected
if (length(break_indices) > 1) {
  # Ensure matching lengths for breaks and labels
  breaks_vec <- c(min(ndvi_df_drought$day), break_indices, max(ndvi_df_drought$day))
  labels_vec <- paste("Segment", seq_len(length(breaks_vec) - 1))
  # Assign segments to the dataset
  ndvi_df_drought$segment <- cut(ndvi_df_drought$day, 
                                 breaks = breaks_vec, 
                                 labels = labels_vec)
  # Fit models for each segment
  lm_segment1 <- lm(ndvi ~ 1, data = ndvi_df_drought[ndvi_df_drought$segment == "Segment 1", ])
  lm_segment2 <- lm(ndvi ~ 1, data = ndvi_df_drought[ndvi_df_drought$segment == "Segment 2", ])
  lm_segment3 <- lm(ndvi ~ 1, data = ndvi_df_drought[ndvi_df_drought$segment == "Segment 3", ])
  # Plot NDVI time series with breakpoints
  plot(ndvi_df_drought$day, ndvi_df_drought$ndvi, type = "l", col = "black", main = "NDVI with Breakpoints")
  abline(v = ndvi_df_drought$day[break_indices], col = "red", lty = 2)  # Add vertical lines for breakpoints
} else {
  print("Only one breakpoint detected.")
}

# ========================================== STRUCCHANGE::fSTATS MODEL BEETLE INFESTATION =========================================

ndvi_df_beetle$ndvi_smoothed <- rollmean(ndvi_df_beetle$ndvi, k = 5, fill = NA)

#fstats_model <- Fstats(ndvi_df_drought$ndvi ~ 1, data = ndvi_df_drought) #This formula specifies a constant mean model (changes in the mean are tested).
fstats_model <- Fstats(ndvi_df_beetle$ndvi_smoothed ~ ndvi_df_beetle$day, data = ndvi_df_beetle) #This formula specifies a constant mean model (changes in the mean are tested).
#fstats_model <- Fstats(ndvi_df_beetle$ndvi ~ 1, data = ndvi_df_beetle) #This formula specifies a constant mean model (changes in the mean are tested).
# Plot the fstats model
plot(fstats_model,main = "F-Statistics for Structural Changes - Drought Event")
# Print the fstats model
print(fstats_model)  # Print the fstats model
# Extract the breakpoints
breakpoints_fstats <- breakpoints(fstats_model)
# Print the breakpoints
print(breakpoints_fstats)  # Print the breakpoints
# Extract the breakpoint observation number
breakpoint_observation <- breakpoints_fstats$breakpoints
breakdates_result <- breakdates(breakpoints_fstats)
# Print the extracted breakpoint
print(breakpoint_observation)
print(breakdates_result)
# Add vertical lines at the breakpoint observation numbers
abline(v = breakdates_result, col = "red", lty = 2, lwd = 2)  # Dashed red lines

# ================================================

# ========================================== BARK BEETLE MODEL =========================================
# Apply the strucchange::fstats function
fstats_model <- Fstats(ndvi_df_beetle$ndvi ~ 1, data = ndvi_df_beetle) #This formula specifies a constant mean model (changes in the mean are tested).
# Plot the fstats model
plot(fstats_model,main = "F-Statistics for Structural Changes - Bark Beetle infestation")
# Print the fstats model
print(fstats_model)  # Print the fstats model
# Extract the breakpoints
breakpoints_fstats <- breakpoints(fstats_model)
# Print the breakpoints
print(breakpoints_fstats)  # Print the breakpoints
# Extract the breakpoint observation number
breakpoint_observation <- breakpoints_fstats$breakpoints
breakdates_result <- breakdates(breakpoints_fstats)
# Print the extracted breakpoint
print(breakpoint_observation)
print(breakdates_result)
# Add vertical lines at the breakpoint observation numbers
abline(v = breakdates_result, col = "red", lty = 2, lwd = 2)  # Dashed red lines


# ======================================== BCP MODEL =========================================

library(bcp)
library(grid)

# Remove NA values from the ndvi data
clean_ndvi <- na.omit(ndvi_df_beetle$ndvi)

bcp_model <- bcp(clean_ndvi)
#bcp_model <- bcp(ndvi_df_beetle$ndvi)
# Plot the bcp model
plot(bcp_model, main="Univariate Change Point - Bark Beettle")
legacyplot(bcp_model)
# Print the bcp model
print(bcp_model)
# Extract posterior probabilities
posterior_probs <- bcp_model$posterior.prob
# Identify the highest probability points (e.g., threshold = 1)
high_prob_indices <- which(posterior_probs > 0.95 & posterior_probs < 0.999)
# Print the indices of highest probability change points
print(high_prob_indices)
# Base plot
plot(bcp_model)
# Plot the bcp model
legacyplot(bcp_model)
# Plot the highest probability points along the x-axis
points(high_prob_indices, rep(0, length(high_prob_indices)), col = "red", pch = 19, cex = 1.5)  # Red dots at y = 0


bcp_model <- bcp(ndvi_df_drought$ndvi)
# Plot the bcp model
plot(bcp_model, main="Univariate Change Point - Drought event")
legacyplot(bcp_model)
# Print the bcp model
print(bcp_model)
# Extract posterior probabilities
posterior_probs <- bcp_model$posterior.prob
# Identify the highest probability points (e.g., threshold = 1)
high_prob_indices <- which(posterior_probs > 0.96)
# Print the indices of highest probability change points
print(high_prob_indices)
# Base plot
plot(bcp_model)

legacyplot(bcp_model)
# Plot the highest probability points along the x-axis
points(high_prob_indices, rep(0, length(high_prob_indices)), col = "red", pch = 19, cex = 1.5)  # Red dots at y = 0


# ========================================= ECP MODEL =========================================
library(ecp)
#library(ggplot2)
#library(ggplotify)
#library(grid)
#library(gridExtra)
#library(ggpubr)


# ndvi_df_beetle ============================================================================
#The ecp package is primarily built around functions like e.divisive() and e.cp3o() for change point analysis. There is no function explicitly called ecp() in the package.

#ndvi_df_drought$ndvi_smoothed <- rollmean(ndvi_df_drought$ndvi, k = 5, fill = NA)

# e.divisive(): For divisive hierarchical estimation of multiple change points.
ecp_model_divisive_beetle <- e.divisive(as.matrix(ndvi_df_beetle$ndvi), min.size = 10, alpha = 0.1)
change_points_beetle <- ecp_model_divisive_beetle$estimates
print(change_points_beetle)  # Print the detected change points
# Create the NDVI plot
plot(ndvi_df_beetle$day, ndvi_df_beetle$ndvi, type = "l", col = "black", main = "NDVI with Change Points - Beetle Infestation")
print(ecp_model_divisive)

# Add vertical lines for each change point
for (cp in change_points_beetle[-1]) {  # Exclude the first point (start of series)
  abline(v = ndvi_df_beetle$day[cp], col = "red", lty = 2)
}

# Quantify Magnitudes of change

# Calculate the magnitude of change between segments
segment_means_beetle <- sapply(seq_along(change_points_beetle[-1]), function(i) {
  start <- change_points_beetle[i]
  end <- change_points_beetle[i + 1] - 1
  mean(ndvi_df_beetle$ndvi[start:end])
})

# Calculate magnitude of changes
magnitude_changes_beetle <- abs(diff(segment_means_beetle))

# Find the most significant change points (e.g., top 2 largest changes)
significant_indices_beetle <- order(magnitude_changes_beetle, decreasing = TRUE)[1:2]
significant_change_points_beetle <- change_points_beetle[significant_indices_beetle + 1]  # Adjust for indices

# Plot NDVI time series
plot(ndvi_df_beetle$day, ndvi_df_beetle$ndvi, type = "l", col = "black", main = "NDVI with Significant Change Points - Beetle Infestation")

# Highlight all change points in red
for (cp in change_points_beetle[-1]) {
  abline(v = ndvi_df_beetle$day[cp], col = "red", lty = 2)
}

# Highlight significant change points in blue
for (scp in significant_change_points_beetle) {
  abline(v = ndvi_df_beetle$day[scp], col = "blue", lty = 1, lwd = 2)  # Solid blue lines
}

legend("topright", legend = c("Change Points", "Significant Points"), col = c("red", "blue"), lty = c(2, 1), lwd = c(1, 2))

# Highlight Posterior probabilities

# Define range for significant probabilities
significant_indices_beetle <- which(posterior_probs > 0.85)
print(significant_indices_beetle)  # Print the indices of significant probabilities
# Highlight these points on the plot
points(ndvi_df_beetle$day[significant_indices_beetle], ndvi_df_beetle$ndvi[significant_indices_beetle], col = "blue", pch = 20, cex = 2)
# Add labels next to each significant point
text(
  x = ndvi_df_beetle$day[significant_indices],        # x-coordinates
  y = ndvi_df_beetle$ndvi[significant_indices],       # y-coordinates
  labels = significant_indices,                       # Labels (indices)
  pos = 4,                                            # Position: to the right of the point
  col = "blue",                                       # Text color
  cex = 0.8                                           # Text size
)

# ============= results do not show significant changes ============================



# ================================ndvi_df_drought ============================================================================
#The ecp package is primarily built around functions like e.divisive() and e.cp3o() for change point analysis. There is no function explicitly called ecp() in the package.

#ndvi_df_drought$ndvi_smoothed <- rollmean(ndvi_df_drought$ndvi, k = 5, fill = NA)

# e.divisive(): For divisive hierarchical estimation of multiple change points.
ecp_model_divisive_drought <- e.divisive(as.matrix(ndvi_df_drought$ndvi), min.size = 10)
change_points_drought <- ecp_model_divisive_drought$estimates
print(change_points_drought)  # Print the detected change points
# Create the NDVI plot
plot(ndvi_df_drought$day, ndvi_df_drought$ndvi, type = "l", col = "black", main = "NDVI with Change Points - Drought Event")

# Add vertical lines for each change point
for (cp in change_points_drought[-1]) {  # Exclude the first point (start of series)
  abline(v = ndvi_df_drought$day[cp], col = "red", lty = 2)
}

# Quantify Magnitudes of change

# Calculate the magnitude of change between segments
segment_means_drought <- sapply(seq_along(change_points_drought[-1]), function(i) {
  start <- change_points_drought[i]
  end <- change_points_drought[i + 1] - 1
  mean(ndvi_df_drought$ndvi[start:end])
})

# Calculate magnitude of changes
magnitude_changes_drought <- abs(diff(segment_means_drought))

# Find the most significant change points (e.g., top 2 largest changes)
significant_indices_drought <- order(magnitude_changes_drought, decreasing = TRUE)[1:2]
significant_change_points_drought <- change_points[significant_indices_drought + 1]  # Adjust for indices

# Plot NDVI time series
plot(ndvi_df_drought$day, ndvi_df_drought$ndvi, type = "l", col = "black", main = "NDVI with Significant Change Points - Drought Event")

# Highlight all change points in red
for (cp in change_points_drought[-1]) {
  abline(v = ndvi_df_drought$day[cp], col = "red", lty = 2)
}

# Highlight significant change points in blue
for (scp in significant_change_points_drought) {
  abline(v = ndvi_df_drought$day[scp], col = "blue", lty = 1, lwd = 2)  # Solid blue lines
}

legend("topright", legend = c("Change Points", "Significant Points"), col = c("red", "blue"), lty = c(2, 1), lwd = c(1, 2))

# Highlight Posterior probabilities

# Define range for significant probabilities
significant_indices_drought <- which(posterior_probs > 0.85)
print(significant_indices_drought)  # Print the indices of significant probabilities
# Highlight these points on the plot
points(ndvi_df_drought$day[significant_indices_drought], ndvi_df_drought$ndvi[significant_indices_drought], col = "blue", pch = 20, cex = 2)
# Add labels next to each significant point
text(
  x = ndvi_df_drought$day[significant_indices],        # x-coordinates
  y = ndvi_df_drought$ndvi[significant_indices],       # y-coordinates
  labels = significant_indices,                       # Labels (indices)
  pos = 4,                                            # Position: to the right of the point
  col = "blue",                                       # Text color
  cex = 0.8                                           # Text size
)



# ========================================= CHANGE POINT MODEL =========================================

library(changepoint)

# Apply the changepoint detection model using the PELT method
# Remove NA values
ndvi_df_beetle <- ndvi_df_beetle[!is.na(ndvi_df_beetle$ndvi), ]
# cpt model: meanvar uses mean and variance simultaneously to estimate changes.
# PELT = Pruned Exact Linear Time. Reduce computational complexity. Usuefull for large datasers
# BIC = Bayesian Information Criterion. Penalizes the number of change points.It prevent overfitting
cpt_model <- cpt.mean(ndvi_df_beetle$ndvi, method = "PELT", penalty = "BIC")

# Print the summary of detected change points
summary(cpt_model)
# Plot the detected change points
plot(cpt_model, main = "Change Point Detection using PELT Method")
# Extract the change points
cpt_points <- cpts(cpt_model)
print(cpt_points)  # Print the detected change points
plot(cpt_model, main = "Change Point Detection using PELT Method")
# Highlight the detected change points
abline(v = cpt_points, col = "red", lty = 2, lwd = 2)  # Red dashed line for change points
legend("topright", legend = "Change Points", col = "red", lty = 2, lwd = 2)


# Extract the NDVI values
ndvi_values <- ndvi_df_beetle$ndvi

# Initialize lists to store segment statistics
mean_differences <- numeric(length(cpt_points) - 1)  # Differences in means between segments

# Calculate the mean difference for each change point
for (i in seq_along(cpt_points)[-1]) {
  segment_1 <- ndvi_values[(cpt_points[i - 1] + 1):cpt_points[i]]
  segment_2 <- ndvi_values[(cpt_points[i] + 1):ifelse(i == length(cpt_points), length(ndvi_values), cpt_points[i + 1])]
  mean_differences[i - 1] <- abs(mean(segment_1) - mean(segment_2))
}

# Find the index of the most significant change point (largest mean difference)
most_significant_index <- which.max(mean_differences)
most_significant_change_point <- cpt_points[most_significant_index]

# Print the most significant change point
print(most_significant_change_point)
#Add to the plot the most significant change point as a vertical line
abline(v = most_significant_change_point, col = "black", lty = 2, lwd = 4)  # Purple dashed line for most significant change point

#========== cpt drought test ===================================================================
# Apply the changepoint detection model using the PELT method
# Remove NA values
ndvi_df_drought <- ndvi_df_drought[!is.na(ndvi_df_drought$ndvi), ]
# cpt model: meanvar uses mean and variance simultaneously to estimate changes.
# PELT = Pruned Exact Linear Time. Reduce computational complexity. Usuefull for large datasers
# BIC = Bayesian Information Criterion. Penalizes the number of change points.It prevent overfitting
cpt_model_drought <- cpt.meanvar(ndvi_df_drought$ndvi, method = "PELT", penalty = "BIC")

# Print the summary of detected change points
summary(cpt_model_drought)
# Plot the detected change points
plot(cpt_model_drought, main = "Change Point Detection using PELT Method")
# Extract the change points
cpt_points <- cpts(cpt_model_drought)
print(cpt_points)  # Print the detected change points
plot(cpt_model_drought, main = "Change Point Detection using PELT Method")
# Highlight the detected change points
#abline(v = cpt_points, col = "red", lty = 2, lwd = 2)  # Red dashed line for change points
legend("topright", legend = "Change Points", col = "red", lty = 2, lwd = 2)


# Initialize lists to store segment statistics
mean_differences <- numeric(length(cpt_points) - 1)  # Differences in means between segments

# Calculate the mean difference for each change point
for (i in seq_along(cpt_points)[-1]) {
  segment_1 <- ndvi_values[(cpt_points[i - 1] + 1):cpt_points[i]]
  segment_2 <- ndvi_values[(cpt_points[i] + 1):ifelse(i == length(cpt_points), length(ndvi_values), cpt_points[i + 1])]
  mean_differences[i - 1] <- abs(mean(segment_1) - mean(segment_2))
}

# Find the index of the most significant change point (largest mean difference)
most_significant_index <- which.max(mean_differences)
most_significant_change_point <- cpt_points[most_significant_index]

# Print the most significant change point
print(most_significant_change_point)
#Add to the plot the most significant change point as a vertical line
abline(v = most_significant_change_point, col = "black", lty = 2, lwd = 4)  # Purple dashed line for most significant change point





# ========================== TSMCP model - it is not working, package has not been updated ==============================
# Load the required library
library(TSMCP)
library(grpreg)
# Apply the tsmcp function to my dataset ndvi_df_beetle$ndvi 
tsmcp_model <- tsmcp(ndvi_df_beetle$ndvi, m = 1, n = 1, alpha = 0.05)
# Print the summary of the tsmcp model
summary(tsmcp_model)
# Plot the tsmcp model
plot(tsmcp_model, main = "Change Point Detection using tsmcp")
# Extract the change points
cpt_points_tsmcp <- tsmcp_model$cpts
print(cpt_points_tsmcp)  # Print the detected change points
# Plot the original NDVI time series
plot(ndvi_df_beetle$ndvi, type = "l", col = "blue", 
     xlab = "Time Index", ylab = "NDVI", 
     main = "NDVI Time Series with Change Points")
# Highlight the detected change points
abline(v = cpt_points_tsmcp, col = "red", lty = 2, lwd = 2)  # Red dashed line for change points
legend("topright", legend = "Change Points", col = "red", lty = 2, lwd = 2)
# Extract the NDVI values
ndvi_values <- ndvi_df_beetle$ndvi
# Initialize lists to store segment statistics
mean_differences <- numeric(length(cpt_points_tsmcp) - 1)  # Differences in means between segments
# Calculate the mean difference for each change point
for (i in seq_along(cpt_points_tsmcp)[-1]) {
  segment_1 <- ndvi_values[(cpt_points_tsmcp[i - 1] + 1):cpt_points_tsmcp[i]]
  segment_2 <- ndvi_values[(cpt_points_tsmcp[i] + 1):ifelse(i == length(cpt_points_tsmcp), length(ndvi_values), cpt_points_tsmcp[i + 1])]
  mean_differences[i - 1] <- abs(mean(segment_1) - mean(segment_2))
}
# Find the index of the most significant change point (largest mean difference)
most_significant_index <- which.max(mean_differences)
most_significant_change_point <- cpt_points_tsmcp[most_significant_index]
# Print the most significant change point
print(most_significant_change_point)
# Add to the plot the most significant change point as a vertical line
abline(v = most_significant_change_point, col = "purple", lty = 2, lwd = 2)  # Purple dashed line for most significant change point


# ================= cpm model ==========================
# cpm - bark beetle infestation

install.packages("cpm")  # Install the cpm package
library(cpm)

# cpm - drough stress =========================================

#cpm_model <- processStream(ndvi_df_drought$ndvi, cpmType = "Student", ARL0 = 1000)
cpm_model <- processStream(ndvi_df_drought$ndvi, cpmType = "Bartlett", ARL0 = 500)
#cpm_model <- processStream(ndvi_df_drought$ndvi, cpmType = "Cramer-von-Mises", ARL0 = 500)
# Print the detected change points
print(cpm_model)
plot(cpm_model$x)  # Plot the cpm model
change_points <- cpm_model$changePoints
print(change_points)  # Display all detected change points

mean_differences <- numeric(length(change_points) - 1)  # Initialize differences

for (i in seq_along(change_points)[-1]) {
  segment_1 <- ndvi_df_drought$ndvi[change_points[i - 1]:change_points[i] - 1]
  segment_2 <- ndvi_df_drought$ndvi[change_points[i]:ifelse(i == length(change_points), length(ndvi_df_drought$ndvi), change_points[i + 1] - 1)]
  mean_differences[i - 1] <- abs(mean(segment_1) - mean(segment_2))
}

# Rank by significance (largest mean differences)
top_two_indices <- order(mean_differences, decreasing = TRUE)[1:1]
most_significant_change_points <- change_points[top_two_indices]

print(most_significant_change_points)  # Display top two most significant change points

# Plot the NDVI data
plot(ndvi_df_drought$ndvi, type = "l", col = "black", main = "NDVI with Most Significant Change Points")
# Add vertical lines for all change points
abline(v = change_points, col = "red", lty = 2)  # All change points
# Add vertical lines for the two most significant change points
abline(v = most_significant_change_points, col = "blue", lty = 1, lwd = 2)  # Highlight most significant ones
# Add the specific vertical green line at x = 1231
abline(v = 833, col = "green", lty = 3, lwd = 4)  # Dashed green line with thicker width
abline(v = 1613, col = "green", lty = 3, lwd = 4)  # Dashed green line with thicker width
# Add text labels next to the most significant change points
for (i in seq_along(most_significant_change_points)) {
  text(x = most_significant_change_points[i],        # X-coordinate (change point)
       y = max(ndvi_df_drought$ndvi, na.rm = TRUE) - (i * 0.05),  # Y-coordinate (adjust height as needed)
       labels = paste("Point:", most_significant_change_points[i]),  # Label text
       pos = 4,                                      # Position (4 = right)
       col = "blue",                                 # Text color
       cex = 0.8)                                    # Text size
}
# Add a legend for clarity
legend("topright", legend = c("All Change Points", "Most Significant"), 
       col = c("red", "blue"), lty = c(2, 1), lwd = c(1, 2))



# cpm - bark beetle infestation =========================================

#cpm_model <- processStream(ndvi_df_drought$ndvi, cpmType = "Student", ARL0 = 1000)
cpm_model <- processStream(ndvi_df_beetle$ndvi, cpmType = "Bartlett", ARL0 = 500)
#cpm_model <- processStream(ndvi_df_drought$ndvi, cpmType = "Cramer-von-Mises", ARL0 = 500)
# Print the detected change points
print(cpm_model)
plot(cpm_model$x)  # Plot the cpm model
change_points <- cpm_model$changePoints
print(change_points)  # Display all detected change points

mean_differences <- numeric(length(change_points) - 1)  # Initialize differences

for (i in seq_along(change_points)[-1]) {
  segment_1 <- ndvi_df_beetle$ndvi[change_points[i - 1]:change_points[i] - 1]
  segment_2 <- ndvi_df_beetle$ndvi[change_points[i]:ifelse(i == length(change_points), length(ndvi_df_drought$ndvi), change_points[i + 1] - 1)]
  mean_differences[i - 1] <- abs(mean(segment_1) - mean(segment_2))
}

# Rank by significance (largest mean differences)
top_two_indices <- order(mean_differences, decreasing = TRUE)[1:2]
most_significant_change_points <- change_points[top_two_indices]

print(most_significant_change_points)  # Display top two most significant change points

# Plot the NDVI data
plot(ndvi_df_beetle$ndvi, type = "l", col = "black", main = "NDVI with Most Significant Change Points")
# Add vertical lines for all change points
abline(v = change_points, col = "red", lty = 2)  # All change points
# Add vertical lines for the two most significant change points
abline(v = most_significant_change_points, col = "blue", lty = 1, lwd = 2)  # Highlight most significant ones
# Add the specific vertical green line at x = 1231
abline(v = 1231, col = "green", lty = 3, lwd = 2)  # Dashed green line with thicker width
# Add text labels next to the most significant change points
for (i in seq_along(most_significant_change_points)) {
  text(x = most_significant_change_points[i],        # X-coordinate (change point)
       y = max(ndvi_df_beetle$ndvi, na.rm = TRUE) - (i * 0.05),  # Y-coordinate (adjust height as needed)
       labels = paste("Point:", most_significant_change_points[i]),  # Label text
       pos = 4,                                      # Position (4 = right)
       col = "blue",                                 # Text color
       cex = 0.8)                                    # Text size
}
# Add a legend for clarity
legend("topright", legend = c("All Change Points", "Most Significant"), 
       col = c("red", "blue"), lty = c(2, 1), lwd = c(1, 2))







































