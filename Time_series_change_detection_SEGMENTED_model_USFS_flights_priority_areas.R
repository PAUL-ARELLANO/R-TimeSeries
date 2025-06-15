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
#setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/Stresses_dataset_28-Apr-25/SUBSET_testing") # Change this path to your directory where NDVI (or other) indexes are stored.
setwd("C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/Synthetic_datasets/PLANET_NDVI/3_0/subset") # This path contains samples extracted from Revanth-Planet weekly time series for specific pixels.

# Retrieve CSV files that start with "NDVI"
#file_list <- list.files(getwd(), pattern = "^NDVI.*\\.csv$", full.names = TRUE)
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
  if (is.null(fit_segmented$psi) || nrow(fit_segmented$psi) == 0) {
    return(data.frame(file = NA, significant_cp = -1, matrix(rep(-1, npsi), nrow = 1, dimnames = list(NULL, paste0("change_point_", 1:npsi)))))
  } else {
    change_points <- fit_segmented$psi[, "Est."]
    
    # Ensure output format for exactly npsi breakpoints, filling missing values with -1
    change_points <- c(change_points, rep(-1, max(0, npsi - length(change_points))))  
    
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
    significant_cp <- ifelse(all(is.na(ndvi_values)), -1, change_points[which.max(ndvi_values)])
    
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
  # Read dataset
  data <- fread(file)
  
  # Ensure NDVI column exists and is non-empty
  if (!("ndvi" %in% colnames(data))) {
    print(paste("Skipping file:", file, "- 'ndvi' column missing"))
    return(NULL)  # Skip this file
  }
  
  if (all(is.na(data$ndvi)) || length(data$ndvi) == 0) {
    print(paste("Skipping file:", file, "- 'ndvi' column is empty"))
    return(NULL)  # Skip this file
  }  
  
  # Fit segmented model and extract change points
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

# Plot histogram of significant change points
hist_plot <- ggplot(results_df, aes(x = significant_cp)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Significant Change Points",
       x = "Significant Change Point", y = "Frequency") +
  theme_minimal()

# Display plot
print(hist_plot)

# Create cumulative distribution plot
cumulative_plot <- ggplot(results_df, aes(x = significant_cp)) +
  stat_ecdf(geom = "step", color = "blue") +
  labs(title = "Cumulative Distribution of Significant Change Points",
       x = "Significant Change Point", y = "Cumulative Probability") +
  theme_minimal()

# Display plot
print(cumulative_plot)


