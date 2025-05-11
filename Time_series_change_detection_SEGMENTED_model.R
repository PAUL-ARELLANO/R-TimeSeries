# ============================================================
# NDVI Change-Point Detection Using Segmented Regression
# ============================================================

# This script processes NDVI time-series data to detect structural shifts
# in vegetation trends using segmented regression. It automates the 
# retrieval, processing, and analysis of multiple CSV files containing NDVI values.

# ============================================================
# 1. Load Required Libraries
# ============================================================

library(segmented)  # For detecting change points using segmented regression
library(data.table)  # Efficient handling of large datasets
library(furrr)  # Parallel processing to speed up execution
library(ggplot2)  # Visualization of detected change points

# ============================================================
# 2. Retrieve NDVI Data Files
# ============================================================

# Set working directory where CSV files are stored
setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets")

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
# 4. Change-Point Detection Using Segmented Regression
# ============================================================

# Function to apply segmented regression and detect structural shifts
fit_segmented_model <- function(data) {
  # Ensure dataset contains required columns
  if (!all(c("ndvi", "day") %in% colnames(data))) {
    stop("Data must contain 'ndvi' and 'day' columns.")
  }
  
  # Apply NDVI preprocessing
  data$ndvi <- preprocess_ndvi(data$ndvi)
  
  # Fit the base linear model
  fit_lm <- lm(ndvi ~ 1 + day, data = data)  # Establish linear trend
  
  # Apply segmented regression to detect change points
  fit_segmented <- tryCatch({
    segmented(fit_lm, seg.Z = ~ day, npsi = 1, 
              control = seg.control(it.max = 50, fix.npsi = TRUE, display = FALSE))
  }, error = function(e) {
    return(NULL)  # Handle segmentation failures
  })
  
  # Ensure valid change-point extraction
  if (is.null(fit_segmented) || is.null(fit_segmented$psi)) {
    return(data.frame(change_point = NA, x_value = NA))  # Return missing values if detection fails
  } else {
    change_point <- fit_segmented$psi[, "Est."]
    x_value <- data$day[which.min(abs(data$day - change_point))]  # Locate closest day to detected change point
    
    return(data.frame(change_point = change_point, x_value = x_value))
  }
}

# ============================================================
# 5. Process Multiple NDVI Datasets in Parallel
# ============================================================

# Function to process each file and extract change points
process_file <- function(file) {
  data <- fread(file)  # Read dataset
  
  # Fit segmented model and extract change points
  change_point_data <- fit_segmented_model(data)
  
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

# Save results to a CSV file (commented out for flexibility)
# write.csv(results_df, "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/segmented_model/results_segmented_model.csv", row.names = FALSE)

# Display results
print(results_df)

# ============================================================
# 6. Visualization of Results
# ============================================================

# Histogram of detected change points
hist_plot_x <- ggplot(results_df, aes(x = x_value)) +
  geom_histogram(binwidth = 5, fill = "red", color = "black", alpha = 0.7) +
  labs(title = "Segmented Model: Histogram of Corresponding X Values",
       x = "X Value (Day)", y = "Frequency") +
  theme_minimal()
print(hist_plot_x)

# Cumulative distribution plot of detected change points
cumulative_plot <- ggplot(results_df, aes(x = change_point)) +
  stat_ecdf(geom = "step", color = "blue") +
  labs(title = "Cumulative Distribution of Change-Point Estimates",
       x = "Change-Point Value", y = "Cumulative Probability") +
  theme_minimal()
print(cumulative_plot)
