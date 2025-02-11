library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggstatsplot)
library(patchwork)

file_path <- "!/R_scripts/final_peak_analysis_summary_all.xlsx"


colnames(data) <- tolower(colnames(data))
data <- data %>%
  rename(
    peak_caller = method,
    histone_mark = histone,
    sample = sample
  )

# expected_columns <- c("peak_caller", "histone_mark", "sample", "total_peaks")


data$peak_caller <- as.factor(data$peak_caller)

histone_marks <- unique(data$histone_mark)

for (mark in histone_marks) {
  mark_data <- data %>% filter(histone_mark == mark)
  plot <- ggbetweenstats(
    data = mark_data,
    x = peak_caller,
    y = total_peaks,
    fill = peak_caller,
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    p.adjust.method = "fdr",
    title = paste("Distribution of Peak Numbers for", mark, "Across Different Peak Callers"),
    xlab = "Peak Caller",
    ylab = "Total Number of Peaks",
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

