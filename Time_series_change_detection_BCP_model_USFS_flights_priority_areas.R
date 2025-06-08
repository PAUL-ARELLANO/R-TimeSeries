# ============================================================
# NDVI Change-Point Detection Using Bayesian Change-Point Model (BCP)
# ============================================================

# This script processes NDVI time-series data to detect structural shifts
# in vegetation trends using Bayesian change-point analysis (`bcp` package).
# It automates the retrieval, processing, and analysis of multiple CSV files.

# ============================================================
# 1. Load Required Libraries
# ============================================================

library(bcp)  # Bayesian Change-Point Detection
library(data.table)  # Efficient handling of large datasets
library(furrr)  # Parallel processing to speed up execution
library(ggplot2)  # Visualization of detected change points
library(viridis)  
library(zoo)
library(dplyr)  # Load the package
library(tidyr)  # Load the package
# ============================================================
# 2. Retrieve NDVI Data Files
# ============================================================

# Set working directory where CSV files are stored
setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25/SUBSET_testing")

# Retrieve CSV files that start with "NDVI"
file_list <- list.files(getwd(), pattern = "^NDVI.*\\.csv$", full.names = TRUE)
print(file_list)

# ============================================================
# 3. Preprocessing NDVI Values
# ============================================================

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

# Function to process each file and extract change points
  process_file <- function(file) {
  data <- fread(file)  # Read dataset
  
  # Fit BCP model and extract detected change points with probabilities
  change_point_data <- fit_bcp_model(data)
  
  # Store results with filename
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
write.csv(results_df, "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/bcp_model/results_bcp_model_Synthetic_datasets_28-Apr-25_folder_multiple_change_Points_threshold_097_099_subset_max_probability.csv", row.names = FALSE)

# Display results
print(results_df)

# ============================================================
# 6. Visualization of Results
# ============================================================

results_df <- read.csv("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/bcp_model/results_bcp_model_Synthetic_datasets_folder_multiple_change_Points.csv", header = TRUE)


# Histogram of detected change points
#hist_plot_x <- ggplot(results_df, aes(x = x_value)) +
#  geom_histogram(binwidth = 10, fill = "red", color = "black", alpha = 0.7) +
#  labs(title = "BCP Model: Histogram of Corresponding X Values",
#       x = "X Value (Day)", y = "Frequency") +
#  theme_minimal()
#print(hist_plot_x)
#
# ==============================
#
# Create histogram
hist_plot <- ggplot(results_df, aes(x = corresponding_change_point)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of the change_points with highest_max_posterior_probability",
       x = "corresponding_change_point", y = "Frequency") +
  theme_minimal()

# Display plot
print(hist_plot)





#
# Select only the first three change point columns
change_point_cols <- results_df %>%
  select(matches("^change_point_[1-3]$"))  # Select columns "change_point_1", "change_point_2", and "change_point_3"
#  select(matches("^change_point_[1]$"))  # Select columns "change_point_1", "change_point_2", and "change_point_3"

# Reshape data into long format for ggplot
long_results_df <- change_point_cols %>%
  pivot_longer(cols = everything(), names_to = "change_point_type", values_to = "x_value")

# Histogram with distinct colors for each change point type
hist_plot_x <- ggplot(long_results_df, aes(x = x_value, fill = change_point_type)) +
  geom_histogram(binwidth = 10, color = "black", alpha = 0.7, position = "identity") +
  scale_fill_viridis_d(option = "plasma") +  # Use Viridis discrete color scale
  labs(title = "BCP Model: Histogram of First Three Change Points",
       x = "Change Point (Day)", y = "Frequency", fill = "Change Point Type") +
  theme_minimal()

print(hist_plot_x)

# Cumulative distribution plot of detected change points
cumulative_plot <- ggplot(results_df, aes(x = change_point)) +
  stat_ecdf(geom = "step", color = "blue") +
  labs(title = "BCP Model: Cumulative Distribution of Change-Point Estimates",
       x = "Change-Point Value", y = "Cumulative Probability") +
  theme_minimal()
print(cumulative_plot)
