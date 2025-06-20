library(ggtern)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)

data <- read.csv("all_performance_metrics.csv")
create_ternary_data <- function(data, metric_col) {
  histone_marks <- c("H3K27ac", "H3K4me3", "H3K27me3")
  
  # Normalize the data for ternary plot
  normalized_data <- data %>%
    filter(HistoneMark %in% histone_marks) %>%
    select(Sample, HistoneMark, Method, all_of(metric_col)) %>%
    group_by(Sample, Method) %>%
    mutate(
      metric_sum = sum(.data[[metric_col]], na.rm = TRUE),
      normalized_value = ifelse(metric_sum == 0, 0, .data[[metric_col]] / metric_sum)
    ) %>%
    ungroup() %>%
    select(Sample, HistoneMark, Method, normalized_value) %>%
    tidyr::pivot_wider(
      names_from = HistoneMark, 
      values_from = normalized_value,
      values_fill = 0
    ) %>%
    mutate(Method = case_when(
      Method == "GOPEAKS" ~ "GoPeaks",
      Method == "LANCEOTRON" ~ "LanceOtron",
      TRUE ~ Method
    ))
  return(normalized_data)
}

create_ternary_plot <- function(data, title) {
  ggtern(data = data, aes(x = H3K27ac, y = H3K4me3, z = H3K27me3)) +
    geom_point(aes(color = Method), size = 5, alpha = 0.8) +
    scale_color_manual(
      name = "Method",
      values = c("GoPeaks" = "#1f78b4", "LanceOtron" = "#ffcc00", 
                 "MACS2" = "#4d4d4d", "SEACR" = "#e31a1c")
    ) +
    labs(
      title = title,
      x = "H3K27ac",
      y = "H3K4me3",
      z = "H3K27me3"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "bottom",
      plot.margin = margin(20, 35, 20, 20) 
    ) +
    theme_showgrid()
}
data_precision <- create_ternary_data(data, "Precision")
data_recall <- create_ternary_data(data, "Recall")
data_f1 <- create_ternary_data(data, "F1")

plot1 <- create_ternary_plot(data_precision, "Precision")
plot2 <- create_ternary_plot(data_recall, "Recall") 
plot3 <- create_ternary_plot(data_f1, "F1 Score")

combined_plot <- plot1 + plot2 + plot3 + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("ternary_plots_combined.pdf", combined_plot, width = 18, height = 7)