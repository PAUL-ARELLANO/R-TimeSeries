# ============================================================
# NDVI Change-Point Detection Using Bayesian Change-Point Model (BCP)
# ============================================================

# This script processes NDVI time-series data to detect structural shifts
# in vegetation trends using Bayesian change-point analysis (`bcp` package).
# It automates the retrieval, processing, and analysis of multiple CSV files.
#
# Updated at Jun 14, 2025
#
# Paul Arellano, Ph.D.
#
# ============================================================
# 1. Load Required Libraries
# ============================================================

library(bcp)  # Bayesian Change-Point Detection
library(data.table)  # Efficient handling of large datasets
library(furrr)  # Parallel processing to speed up execution
library(ggplot2)  # Visualization of detected change points
library(viridis)  # For color palettes in ggplot2
library(zoo) # For handling time series data where observations do not occur at regular intervals
library(dplyr)  # This package provides a set of functions for data manipulation
library(tidyr)  # This package provides functions for reshaping data
library(readxl)   # For reading Excel files
library(writexl)  # For writing updated Excel files
library(stringr)  # For better string handling
# ============================================================
# 2. Retrieve NDVI Data Files
# ============================================================

# Set working directory where CSV files are stored
#setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25/SUBSET_testing") # Change this path to your directory where NDVI (or other) indexes are stored.
#setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/PLANET_NDVI/edited") # This path contains samples extracted from Revanth-Planet weekly time series for specific pixels.
setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/PLANET_NDVI/3_0/subset") # This path contains samples extracted from Revanth-Planet weekly time series for specific pixels.
# Retrieve CSV files that start with "NDVI"
file_list <- list.files(getwd(), pattern = "\\.xlsx$", full.names = TRUE)
print(file_list)
# ============================================================
# 2.1. Extract "year" and "date" from the "File" column of each record
#
# Function to extract "year" and "date" from the "File" column of each record
extract_year_date <- function(file_column) {
  matches <- str_match(file_column, ".*/([0-9]+)/([0-9]+)/[^/]+$")
  if (!is.na(matches[1])) {
    year <- as.numeric(matches[2])  # Extract the year (6th segment)
    date <- as.numeric(matches[3])  # Extract the date (7th segment)
    return(c(year, date))
  } else {
    return(c(NA, NA))  # Handle missing matches
  }
}

# Process each file individually and save as CSV
for (file in file_list) {
  df <- read_excel(file)  # Read Excel file
  
  # Ensure 'File' and 'NDVI' columns exist
  if (!("File" %in% names(df)) || !("NDVI" %in% names(df))) {
    print(paste("Skipping file:", file, "- Missing required columns"))
    next
  }
  
  # Apply extraction function to each row in the 'File' column
  extracted_values <- t(sapply(df$File, extract_year_date))
  df$year <- extracted_values[,1]
  df$date <- extracted_values[,2]
  
  # Rename columns
  colnames(df)[colnames(df) == "NDVI"] <- "ndvi"
  colnames(df)[colnames(df) == "date"] <- "day"
  
  # Keep only required columns: File, ndvi, year, and date
  df <- df[, c("File", "ndvi", "year", "day")]
  
  # Sort by year and date
  df <- df[order(df$year, df$day), ]  
  
  # Save updated file as CSV in the same directory
  output_file <- file.path(dirname(file), paste0(tools::file_path_sans_ext(basename(file)), "_updated.csv"))
  write.csv(df, output_file, row.names = FALSE)
  
  print(paste("Saved:", output_file))
}

print("Processing complete. Updated CSV files saved in the same folder.")

# ============================================================
# 3. Preprocessing NDVI Values
# ============================================================
#
# Retrieve CSV files generated in previous step
file_list <- list.files(getwd(), pattern = "\\.csv$", full.names = TRUE)
print(file_list)
#

# Function to clean and interpolate NDVI values
preprocess_ndvi <- function(ndvi_series) {
  ndvi_series <- replace(ndvi_series, ndvi_series == 0, NA)  # Convert zero values to NA
  ndvi_series <- na.approx(ndvi_series, rule = 2)  # Interpolate missing values
  return(ndvi_series)
}

# ============================================================
# 4. Change-Point Detection Using Bayesian Change-Point Model
# ============================================================

# Function to apply Bayesian Change-Point model to detect structural shifts
fit_bcp_model <- function(data) {
  # Ensure dataset contains required columns
  if (!all(c("ndvi", "day") %in% colnames(data))) {
    stop("Data must contain 'ndvi' and 'day' columns.")
  }
  
  # Apply NDVI preprocessing
  data$ndvi <- preprocess_ndvi(data$ndvi)
  
  # Apply Bayesian Change-Point Model
  bcp_result <- bcp(data$ndvi, burnin = 1000, mcmc = 2000)  # burnin = 1000: Initial samples discarded to ensure stable estimation. mcmc = 5000: Number of iterations in the Markov Chain Monte Carlo (MCMC) process—higher values improve accuracy
  
  # Extract all change points where posterior probability falls within range
  threshold_min <- 0.97
  threshold_max <- 1
  high_prob_indices <- which(bcp_result$posterior.prob > threshold_min & bcp_result$posterior.prob < threshold_max)
  change_points <- data$day[high_prob_indices]
  change_point_probs <- bcp_result$posterior.prob[high_prob_indices]
 
  # Ensure valid extraction and format results correctly
  if (length(change_points) == 0) {
    return(data.frame(change_point_1 = -1, max_posterior_probability_1 = -1))
  } else {
    strongest_index <- which.max(change_point_probs)
    highest_max_posterior_probability <- change_point_probs[strongest_index]
    corresponding_change_point <- change_points[strongest_index]
  }
  
  # Create a structured dataframe storing change points and probabilities
  max_length <- max(length(change_points), length(change_point_probs))
  change_point_df <- data.frame(matrix(ncol = max_length * 2, nrow = 1))
  
  # Assign column names dynamically
  colnames(change_point_df) <- c(paste0("change_point_", seq_len(max_length)), 
                                 paste0("max_posterior_probability_", seq_len(max_length)))
  
  # Populate change points and probabilities
  change_point_df[1, seq_len(length(change_points))] <- change_points
  change_point_df[1, (max_length + 1):(max_length + length(change_point_probs))] <- change_point_probs
  
  # Append highest max posterior probability and corresponding change point as additional columns
  change_point_df$highest_max_posterior_probability <- highest_max_posterior_probability
  change_point_df$corresponding_change_point <- corresponding_change_point
  
  return(change_point_df)
}

# ============================================================
# 5. Process Multiple NDVI Datasets in Parallel
# ============================================================

process_file <- function(file) {
  data <- fread(file)  # Read dataset
  
  # Ensure NDVI column exists and is non-empty
  if (!("ndvi" %in% colnames(data))) {
    print(paste("Skipping file:", file, "- 'ndvi' column missing"))
    return(NULL)  # Skip this file
  }
  
  if (all(is.na(data$ndvi)) || length(data$ndvi) == 0) {
    print(paste("Skipping file:", file, "- 'ndvi' column is empty"))
    return(NULL)  # Skip this file
  }
  
  # Proceed with Bayesian change-point detection
  change_point_data <- fit_bcp_model(data)
  
  result <- cbind(file = basename(file), change_point_data)
  return(result)
}

# Apply parallel processing for efficiency
plan(multisession, workers = 4)  # Adjust based on system resources

# Run parallelized change-point detection across all files
results <- future_map(file_list, process_file, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Compile results into a single data frame
results_df <- bind_rows(results)

# Ensure correct column order: File name, then the new fields, then the remaining data
column_order <- c("file", "highest_max_posterior_probability", "corresponding_change_point",
                  setdiff(colnames(results_df), c("file", "highest_max_posterior_probability", "corresponding_change_point")))

# Reorder dataset
results_df <- results_df[, column_order]


# Save results to a CSV file
write.csv(results_df, "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/bcp_model/results_bcp_model_Synthetic_datasets_28-Apr-25_folder_multiple_change_Points_threshold_097_099_subset_max_probability.csv", row.names = FALSE)# Change to the your folder

# Display results
print(results_df)

# ============================================================
# 6. Visualization of Results - Histograms and Cumulative Distribution
# ============================================================

results_df <- read.csv("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/bcp_model/results_bcp_model_Synthetic_datasets_folder_multiple_change_Points.csv", header = TRUE)

# Create histogram of detected change points
hist_plot <- ggplot(results_df, aes(x = corresponding_change_point)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of the change_points with highest_max_posterior_probability",
       x = "corresponding_change_point", y = "Frequency") +
  theme_minimal()

# Display plot
print(hist_plot)

# Create a cumulative plot of detected change points
cumulative_plot <- ggplot(results_df, aes(x = corresponding_change_point)) +
  stat_ecdf(geom = "step", color = "blue") +
  labs(title = "BCP Model: Cumulative Distribution of Change-Point Estimates",
       x = "Change-Point Value", y = "Cumulative Probability") +
  theme_minimal()
print(cumulative_plot)
