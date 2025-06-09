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
#setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets")
setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25/SUBSET_testing") # Change this path to your directory where NDVI (or other) indexes are stored.


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
fit_segmented_model <- function(data, npsi = 4) {
  if (!all(c("ndvi", "day") %in% colnames(data))) {
    stop("Data must contain 'ndvi' and 'day' columns.")
  }
  
  # Apply NDVI preprocessing
  data$ndvi <- preprocess_ndvi(data$ndvi)
  
  # Fit base linear model
  fit_lm <- lm(ndvi ~ day, data = data)
  
  # Apply segmented regression
  fit_segmented <- tryCatch({
    segmented(fit_lm, seg.Z = ~ day, npsi = npsi, 
              control = seg.control(it.max = 50, fix.npsi = TRUE, display = FALSE))
  }, error = function(e) {
    return(data.frame(file = NA, significant_cp = NA, matrix(rep(NA, npsi), nrow = 1, dimnames = list(NULL, paste0("change_point_", 1:npsi)))))  
  })
  
  # Extract breakpoints
  if (is.null(fit_segmented$psi)) {
    return(data.frame(file = NA, significant_cp = NA, matrix(rep(NA, npsi), nrow = 1, dimnames = list(NULL, paste0("change_point_", 1:npsi)))))
  } else {
    change_points <- fit_segmented$psi[, "Est."]
    
    # Ensure output format for exactly npsi breakpoints, filling missing values with NA
    change_points <- c(change_points, rep(NA, max(0, npsi - length(change_points))))  
    
    # Identify NDVI values before and after each breakpoint
    ndvi_values <- sapply(change_points, function(cp) {
      idx <- which.min(abs(data$day - cp))  # Find closest index to the change point
      if (idx > 1 && idx < nrow(data)) {
        return(abs(data$ndvi[idx] - data$ndvi[idx - 1]))  # Compute NDVI shift
      } else {
        return(NA)  # Avoid indexing errors
      }
    })
    
    # Determine most significant change point
    significant_cp <- change_points[which.max(ndvi_values)]
    
    return(data.frame(file = NA, significant_cp = significant_cp, matrix(change_points, nrow = 1, dimnames = list(NULL, paste0("change_point_", 1:npsi)))))
  }
}

# ============================================================
# 5. Process Multiple NDVI Datasets in Parallel
# ============================================================

# Function to process files
process_file <- function(file, npsi = 4) {
  file <- as.character(file)  # Ensure proper format
  if (!nzchar(file)) return(NULL)  # Handle empty strings
  
  data <- fread(file)
  change_point_data <- fit_segmented_model(data, npsi)
  
  return(data.frame(file = basename(file), change_point_data))
}

# Run parallelized detection
results <- future_map(file_list, ~process_file(as.character(.x), npsi = 4), .progress = TRUE, .options = furrr_options(seed = TRUE))
results_df <- bind_rows(results)

# Save results
write.csv(results_df, "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/segmented_model/results_segmented_model_USFS_flights_priority_areas.csv", row.names = FALSE)

# Display results
print(results_df)
# ============================================================
# 6. Visualization of Results
# ============================================================

plot_significant_cp <- function(data, significant_cp) {
  # Base NDVI time series plot
  plot(data$day, data$ndvi, type = "l", col = "blue", lwd = 2, 
       xlab = "Day", ylab = "NDVI", main = "Most Significant Change Point")
  
  # Highlight significant breakpoint
  if (!is.na(significant_cp)) {
    abline(v = significant_cp, col = "red", lwd = 2, lty = 2)  # Dashed red line for significant CP
    legend("topright", legend = paste("Significant CP:", significant_cp), col = "red", lty = 2, lwd = 2)
  } else {
    text(mean(data$day), max(data$ndvi), "No Significant Change Point Detected", col = "red", cex = 1.2)
  }
}
# Example usage for one dataset
data_example <- fread("NDVI_noise0.05_trend-0.2_decline0.4_NA0.1_amplitude0.6_phase-0.02.csv")  # Load a specific dataset
plot_significant_cp(data_example, results_df$significant_cp[1])  # Plot first dataset's significant CP

