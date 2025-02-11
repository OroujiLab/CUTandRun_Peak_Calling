# Benchmarking Peak Calling Method for CUT&RUN
# Author: Elias Orouji
# Date: November 4, 2024
# Manuscript Title: "Benchmarking Peak Calling Methods for CUT&RUN"
# Figure Number: [Figure 1b]

# Description: This script generates the plots for the distribution of peak length in three histone marks using four peak calling methods.

# Notes:
# - This script is intended for use with data from mouse brain samples with three different histone marks.
# - The peak calling methods compared include MACS2, SEACR, GoPeaks, and LanceOtron.
# - Ensure the appropriate input files are loaded for each figure, and customize the parameters as needed.


# Set working directory
setwd("~/Desktop/R_scripts/")

# Clean up environment
rm(list = ls())

# Install required packages if not already installed
required_packages <- c("readxl", "dplyr", "ggplot2", "tidyr", "ggstatsplot", "patchwork")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggstatsplot)
library(patchwork)

# Load the data from the Excel file
file_path <- "/Users/elias.orouji/Desktop/R_scripts/final_peak_analysis_summary_all.xlsx"
if (file.exists(file_path)) {
  data <- read_excel(file_path)
} else {
  stop("Error: File does not exist. Please check the file path.")
}

# Rename columns to match the expected format
colnames(data) <- tolower(colnames(data))
data <- data %>%
  rename(
    peak_caller = method,
    histone_mark = histone,
    sample = sample
  )

# Ensure required columns are present
expected_columns <- c("peak_caller", "histone_mark", "sample", "mean_width")
if (!all(expected_columns %in% colnames(data))) {
  stop("Error: Required columns are missing from the data. Please check the data format.")
}

# Convert peak_caller to a factor
data$peak_caller <- as.factor(data$peak_caller)

# Generate and display violin plots for each histone mark separately
histone_marks <- unique(data$histone_mark)

for (mark in histone_marks) {
  mark_data <- data %>% filter(histone_mark == mark)
  plot <- ggbetweenstats(
    data = mark_data,
    x = peak_caller,
    y = mean_width,
    fill = peak_caller,
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    p.adjust.method = "fdr",
    title = paste("Distribution of Peak Width for", mark, "Across Different Peak Callers"),
    xlab = "Peak Caller",
    ylab = "Mean Peak Width",
    ggtheme = ggplot2::theme_minimal(),
    package = "RColorBrewer",
    palette = "Set3",
    type = "parametric",
    plotgrid.args = list(ncol = 2),
    results.subtitle = TRUE
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  print(plot)
}
