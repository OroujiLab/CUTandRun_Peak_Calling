library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggstatsplot)
library(patchwork)

# Set file path and read data
file_path <- "~/Downloads/final_peak_analysis_summary_all.csv"
if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}

# Read the data
data <- read.csv(file_path)
colnames(data) <- tolower(colnames(data))

# Convert method to factor
data$method <- as.factor(data$method)
histone_marks <- unique(data$histone)

# Open PDF device
pdf("figure1_adjusted.pdf", width=11, height=8)

for (mark in histone_marks) {
  mark_data <- data %>% filter(histone == mark)
  
  if (nrow(mark_data) == 0) {
    warning(paste("No data found for histone mark:", mark))
    next
  }
  
  plot <- ggbetweenstats(
    data = mark_data,
    x = method,
    y = total_peaks,
    fill = method,
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    p.adjust.method = "fdr",
    title = paste("Distribution of Peak Numbers for", mark, "Across Different Peak Callers"),
    xlab = "Peak Caller",
    ylab = "Total Number of Peaks",
    ggtheme = ggplot2::theme_minimal(),
    palette = "Set3",
    type = "parametric",
    plotgrid.args = list(ncol = 2),
    results.subtitle = TRUE
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.subtitle = element_text(face = "bold"),
      plot.caption = element_text(face = "bold")
    )
  
  print(plot)
  gc()
}

# Close PDF device
dev.off()