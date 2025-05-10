# ============================================================
# R Script for NDVI Change-Point Detection Using MCP Model
# ============================================================
# Author: Paul Arellano, Ph.D.
# Date: April 26, 2025
# Description:
# This script processes synthetic NDVI datasets and applies a Bayesian 
# change-point model (MCP) to detect structural shifts in vegetation index trends.
# It performs parallel processing for efficiency, preprocesses NDVI values,
# fits the MCP model, extracts change-point estimates, and visualizes results.
#
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

library(dplyr)
library(tidyr)
library(progressr)
library(furrr)
library(mcp)
library(future)
library(future.apply)

library(ggplot2)
library(purrr)
library(ggpubr)
#library(ggplotify)
library(gridExtra)
library(grid)
library(zoo)  # Load the package
library(data.table)

# Set the working directory to the location of your CSV files
setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25")
file_list <- list.files("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25", 
                        pattern = "*.csv", full.names = TRUE)
print(file_list)

# Define the MCP model. GAUSSIAN. Change point model with three segments 
model <- list(
  ndvi ~ 1,  # Initial intercept
  ~ day      # Linear trend before the change point
)

# Function to preprocess NDVI values
preprocess_ndvi <- function(ndvi_series) {
  ndvi_series <- replace(ndvi_series, ndvi_series == 0, NA)  # Convert zero values to NA
  ndvi_series <- na.approx(ndvi_series, rule = 2)  # Linear interpolation
  return(ndvi_series)
}

# Plan the use of multiple cores for parallel processing
plan(multisession)  # Use multiple cores

# Enable progress bar
handlers(global = TRUE)
handlers("txtprogressbar")  # Use a text-based progress bar

# Define the processing function
process_file <- function(file) {
  # Create a progress signal
  progressor <- progressr::progressor(along = file_list)
  progressor()  # Update progress for each file
  
  # Read the current CSV file
  file_path <- file.path(file) # Ensure absolute path is used
#  ndvi_df <- read.csv(file_path, header = TRUE)
  ndvi_df <- fread(file_path, header = TRUE)  # Use fread for faster reading
    print(file_path)  # Debugging step
  
  
  
  # Check for required columns
  if (!all(c("ndvi", "day") %in% colnames(ndvi_df))) {
    warning(paste("File", file, "does not contain required columns: 'ndvi' and 'day'. Skipping."))
    return(data.frame(file_name = file, cp_1_mean = NA))
  }
  
  # Preprocess NDVI values
  ndvi_df$ndvi <- preprocess_ndvi(ndvi_df$ndvi)
  
  # Fit the MCP model to the current dataset
  fit <- mcp(model, data = ndvi_df)
  
  # Extract the summary of the fitted model
  fit_summary <- summary(fit)
  
  # Access the mean of cp_1 and handle missing values
  cp_1_mean <- if ("cp_1" %in% fit_summary$name) {
    fit_summary[fit_summary$name == "cp_1", "mean"]
  } else {
    NA
  }
  
  # Return the results as a dataframe row
  return(data.frame(file_name = file, cp_1_mean = cp_1_mean))
}

# Use progressr with future_lapply to display progress
with_progress({
  results <- future_lapply(file_list, process_file)
})

# Combine the results into a single dataframe
results_df <- do.call(rbind, results)


# Print the results
print(results_df)

# Optionally, save the results to a CSV file
write.csv(results_df, "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/R-Synthetics_TimeSeries_generator/results_MCP_Bark_beetle_summary_cleanedNA0.csv", row.names = FALSE)

# Load the results CSV file
results_df <- read.csv("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/R-Synthetics_TimeSeries_generator/results_MCP_Bark_beetle_summary_cleanedNA0.csv", header = TRUE)

# Create the histogram
hist(results_df$cp_1_mean, 
     main = "MCP model: Histogram of detections", 
     xlab = "cp1_values", 
     col = "lightblue", 
     border = "black", 
     breaks = 20)  # Adjust number of breaks if needed

# Add vertical lines for the confidence interval
#abline(v = mean_cp1, col = "red", lwd = 2)  # Mean line
#abline(v = lower_bound, col = "green", lwd = 2, lty = 2)  # Lower boundary
#abline(v = upper_bound, col = "green", lwd = 2, lty = 2)  # Upper boundary
# Add a vertical line for the BreakPointValue (e.g., 1231)
abline(v = 1231, col = "purple", lwd = 2, lty = 3)  # Dashed purple vertical line
# Define the BreakPointValue and offsets
BreakPointValue <- 1231
offset_plus <- BreakPointValue + 50
offset_minus <- BreakPointValue - 10
abline(v = offset_plus, col = "green", lwd = 2, lty = 2)             # +30
abline(v = offset_minus, col = "green", lwd = 2, lty = 2)            # -30

# Count values within the range
count_in_range <- sum(results_df$cp_1_mean >= offset_minus & results_df$cp_1_mean <= offset_plus, na.rm = TRUE)

# Compute the percentage relative to total samples
total_samples <- length(results_df$cp_1_mean)
percentage_in_range <- (count_in_range / total_samples) * 100

# Print results
print(paste("Number of samples within", offset_minus, "and", offset_plus, ":", count_in_range))
print(paste("Percentage of samples within the range:", round(percentage_in_range, 2), "%"))


# Calculate distance for every record
#
#
#
results_df$Difference <- (results_df$cp_1_mean - 1231)
# Plot bar plot of the differences

#
# remove outliers Difference
results_df <- results_df %>%
  filter(Difference >= -100 & Difference <= 100)

# Add the vertical line for the BreakPointValue = 0
ggplot(results_df, aes(x = Difference)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Histogram of Differences from BreakPointValue",
       x = "Difference from BreakPointValue",
       y = "Frequency") +
  theme_minimal()

# Add cumulative distribution function
ggplot(results_df, aes(x = Difference)) +
  stat_ecdf(geom = "step", color = "blue") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Cumulative Distribution Function of Differences from BreakPointValue",
       x = "Difference from BreakPointValue",
       y = "Cumulative Probability") +
  theme_minimal()
#
# Script end ===========================================================
