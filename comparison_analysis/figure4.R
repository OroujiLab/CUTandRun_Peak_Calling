library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(RColorBrewer)

data <- read.delim("./Transfers/snr_summary_all.tsv", header=TRUE, stringsAsFactors=FALSE)
data2 <- data %>%
  mutate(
    Clean_Sample = case_when(
      grepl("_H3K\\w+3_", Sample) ~ gsub("^.*?_H3K\\w+3_", "", Sample),
      grepl("_H3K27ac_", Sample) ~ gsub("^.*?_H3K27ac_", "", Sample),
      TRUE ~ Sample  
    ),
    Legend_Sample = gsub("_R\\d+$", "", Clean_Sample)
  )

snr_data2 <- data2 %>%
  group_by(Mark, Clean_Sample, Caller) %>%
  summarise(
    peak_mean = mean(Value[Metric == "cnr_peak_Mean"], na.rm = TRUE),
    background_mean = mean(Value[Metric == "cnr_background_Mean"], na.rm = TRUE),
    SD = sd(Value[Metric == "cnr_peak_Mean"], na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Either set a minimum background value
    background_mean = ifelse(background_mean == 0, min(background_mean[background_mean > 0]), background_mean),
    # Or you could set SNR to NA for these cases
    SNR = peak_mean / background_mean,
    SE = ifelse(Count > 1, SD / sqrt(Count), NA),
    Value_Label = round(SNR, 2),
    Legend_Sample = gsub("_R\\d+$", "", Clean_Sample)
  ) %>%
  filter(!is.infinite(SNR))


#snr_data2 <- snr_data2 %>%
#  filter(
#    !Legend_Sample %in% c("H1", "H1_Definitive_Endoderm") &
#      (!Legend_Sample == "K562" | (Legend_Sample == "K562" & Clean_Sample %in% c("K562_R1", "K562_R2")))
#  )

n_samples <- length(unique(snr_data2$Clean_Sample))
colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_samples)

p2 <- ggplot(snr_data2, aes(x = Clean_Sample, y = SNR, fill = Legend_Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.7) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "orange")) +  
  facet_grid(Mark ~ Caller, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Signal-to-Noise Ratio Comparison Across Methods and Samples",
    x = "Sample",
    y = "Signal-to-Noise Ratio",
    fill = "Sample"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 40, r = 10, b = 60, l = 10, unit = "pt"))

ggsave(
  filename = "snr_dataset_filtered2.jpeg",
  plot = p2,
  width = 17,
  height = 10
)

####### for log10

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(RColorBrewer)

data <- read.delim("./Transfers/snr_summary_all.tsv", header=TRUE, stringsAsFactors=FALSE)
data2 <- data %>%
  mutate(
    Clean_Sample = case_when(
      grepl("_H3K\\w+3_", Sample) ~ gsub("^.*?_H3K\\w+3_", "", Sample),
      grepl("_H3K27ac_", Sample) ~ gsub("^.*?_H3K27ac_", "", Sample),
      TRUE ~ Sample  
    ),
    Legend_Sample = gsub("_R\\d+$", "", Clean_Sample)
  )

snr_data2 <- data2 %>%
  group_by(Mark, Clean_Sample, Caller) %>%
  summarise(
    peak_mean = mean(Value[Metric == "cnr_peak_Mean"], na.rm = TRUE),
    background_mean = mean(Value[Metric == "cnr_background_Mean"], na.rm = TRUE),
    SD = sd(Value[Metric == "cnr_peak_Mean"], na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    background_mean = ifelse(background_mean == 0, min(background_mean[background_mean > 0]), background_mean),
    SNR = peak_mean / background_mean,
    SNR_log = log2(SNR),  # log2 transformation
    SE = ifelse(Count > 1, SD / sqrt(Count), NA),
    Value_Label = round(SNR_log, 2),  # Changed to show log2 values
    Legend_Sample = gsub("_R\\d+$", "", Clean_Sample)
  )


#snr_data2 <- snr_data2 %>%
#  filter(
#    !Legend_Sample %in% c("H1", "H1_Definitive_Endoderm") &
#      (!Legend_Sample == "K562" | (Legend_Sample == "K562" & Clean_Sample %in% c("K562_R1", "K562_R2")))
#  )

n_samples <- length(unique(snr_data2$Legend_Sample))
colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_samples)

p2 <- ggplot(snr_data2, aes(x = Clean_Sample, y = SNR_log, fill = Legend_Sample)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.7) +
  facet_grid(Mark ~ Caller, scales = "free_x") +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "orange")) +  
  theme_minimal(base_size = 14) +
  labs(
    title = "Signal-to-Noise Ratio Comparison Across Methods and Samples (log2 scale)",  
    x = "Sample",
    y = "Signal-to-Noise Ratio (log2)",
    fill = "Sample"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 40, r = 10, b = 60, l = 10, unit = "pt"))

ggsave(
  filename = "snr_dataset_filtered_log2.jpeg",
  plot = p2,
  width = 15,
  height = 10
)

write.csv(snr_data2, file = "snr_filtered_plot_data_with_log2.csv", row.names = FALSE)
