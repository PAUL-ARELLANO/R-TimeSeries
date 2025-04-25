# R-TimeSeries
###############################################################################
# Title:       Synthetic NDVI Time Series Generator with Bark Beetle Simulation
# Description: This script generates synthetic NDVI time series data for multiple 
#              scenarios using variable noise, trend, amplitude, seasonal phase 
#              shifts, and a simulated bark beetle infestation event. 
#              It also introduces random NA values and ensures minimum NDVI thresholds.
#
# Author:      Paul Arellano
# Date:        2025-04-25
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
