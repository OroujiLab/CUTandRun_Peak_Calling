library(GenomicRanges)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork) 

file_paths_chip <- list(
  LANCEOTRON = Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/ChiP/VISUALIZATION/results/LANCEOTRON/beds/*H3K*_R*.bed"),
  MACS2 = c(
    Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/ChiP/VISUALIZATION/results/MACS2/*H3K27me3*_R*.broadPeak"),
    Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/ChiP/VISUALIZATION/results/MACS2/*H3K*_R*.narrowPeak")
  ),
  GOPEAKS = Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/ChiP/VISUALIZATION/results/GOPEAKS/*H3K*_R*_peaks.bed"),
  SEACR = Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/ChiP/VISUALIZATION/results/SEACR/*H3K*_R*_stringent*.bed")
)

file_paths_cnr <- list(
  LANCEOTRON = Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/LANCEOTRON/*H3K*_R*.bed"),
  MACS2 = c(
    Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/MACS2/*H3K27me3*_R*.broadPeak"),
    Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/MACS2/*H3K*_R*.narrowPeak")
  ),
  GOPEAKS = Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/GOPEAKS/*H3K*_R*_peaks.bed"),
  SEACR = Sys.glob("/Volumes/Extreme Pro/CUT&RUN/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/SEACR/*H3K*_R*_stringent*.bed")
)

pdf("ChIP_vs_CNR_comparison.pdf", width = 22, height = 14)

count_peaks <- function(file) {
  peaks <- read_delim(file, delim = "\t", col_names = FALSE)
  return(nrow(peaks))
}

standardize_sample_name <- function(sample) {
  # Replace "H1-hESC" with "H1_ESC"
  if (sample == "H1-hESC") {
    return("H1_ESC")
  }
  return(sample)
}

process_files <- function(file_paths, data_type) {
  peak_counts <- data.frame(Method = character(), Sample = character(), 
                            Histone = character(), Replicate = character(), 
                            Peaks = integer(), DataType = character())
  
  for (method in names(file_paths)) {
    for (file in file_paths[[method]]) {
      # Extract file information
      file_name <- basename(file)
      parts <- strsplit(file_name, "_")[[1]]
      histone <- parts[grep("H3K", parts)][1]
      replicate <- parts[grep("^R\\d+", parts)][1]
      histone_index <- which(parts == histone)
      replicate_index <- which(parts == replicate)
      if (any(grepl("H1-hESC", file_name))) {
        sample <- "H1_ESC"  # Standardize to H1_ESC
      } else {
        sample <- paste(parts[(histone_index + 1):(replicate_index - 1)], collapse = "_")
        sample <- standardize_sample_name(sample)  
      }
      num_peaks <- count_peaks(file)
      peak_counts <- rbind(peak_counts, data.frame(Method = method, 
                                                   Sample = sample, 
                                                   Histone = histone, 
                                                   Replicate = replicate, 
                                                   Peaks = num_peaks,
                                                   DataType = data_type))
    }
  }
  
  return(peak_counts)
}

chip_peak_counts <- process_files(file_paths_chip, "ChIP-seq")
cnr_peak_counts <- process_files(file_paths_cnr, "CUT&RUN")

combined_peak_counts <- rbind(chip_peak_counts, cnr_peak_counts)
combined_peak_counts <- combined_peak_counts %>%
  mutate(Histone_Sample = paste(Histone, Sample, sep = "_")) %>%
  mutate(Method_Replicate = paste(Method, Replicate, sep = "_")) %>%
  mutate(Sample_Replicate = paste(Sample, Replicate, sep = "_"),
         Sample_Only = Sample)

combined_peak_counts_filtered <- combined_peak_counts %>%
  filter(
    !Sample %in% c("H1", "H1_Definitive_Endoderm") &
      (!Sample == "K562" | (Sample == "K562" & Replicate %in% c("R1", "R2")))
  )

# Find samples that have both ChIP-seq and CUT&RUN data
samples_in_both <- combined_peak_counts_filtered %>%
  group_by(Method, Histone, Sample, Replicate) %>%
  summarise(DataTypes = n_distinct(DataType), .groups = 'drop') %>%
  filter(DataTypes == 2) %>%
  mutate(Sample_Method_Histone_Rep = paste(Sample, Method, Histone, Replicate, sep = "_"))

# Filter to keep only samples that have both data types
combined_peak_counts_filtered <- combined_peak_counts_filtered %>%
  mutate(Sample_Method_Histone_Rep = paste(Sample, Method, Histone, Replicate, sep = "_")) %>%
  filter(Sample_Method_Histone_Rep %in% samples_in_both$Sample_Method_Histone_Rep)

# Check for and merge duplicates by summing peaks for duplicate entries
combined_peak_counts_merged <- combined_peak_counts_filtered %>%
  group_by(Method_Replicate, Method, Sample_Replicate, Sample_Only, Histone_Only = Histone, DataType) %>%
  summarise(Peaks = sum(Peaks), .groups = 'drop')

# Create a combined method-datatype variable for faceting
combined_peak_counts_merged <- combined_peak_counts_merged %>%
  mutate(Method_DataType = paste(Method, DataType, sep = "_"))

peak_count_plot <- ggplot(combined_peak_counts_merged, aes(x = Sample_Replicate, y = Peaks, fill = DataType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.7) +
  facet_grid(Histone_Only ~ Method, scales = "free_x") +
  scale_fill_manual(values = c("ChIP-seq" = "#4C72B0", "CUT&RUN" = "#DD8452")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Peak Counts: ChIP-seq vs CUT&RUN",
    x = "Sample_Replicate",
    y = "Number of Peaks (log10)",
    fill = "Data Type"
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
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))

print(peak_count_plot)

# Similarly update the calculate_peak_widths function
calculate_peak_widths <- function(file_paths, data_type) {
  peak_widths <- data.frame(Method = character(), Sample = character(), 
                            Histone = character(), Replicate = character(), 
                            Width = numeric(), DataType = character())
  
  for (method in names(file_paths)) {
    for (file in file_paths[[method]]) {
      # Extract file information
      file_name <- basename(file)
      parts <- strsplit(file_name, "_")[[1]]
      histone <- parts[grep("H3K", parts)][1]
      replicate <- parts[grep("^R\\d+", parts)][1]
      histone_index <- which(parts == histone)
      replicate_index <- which(parts == replicate)
      if (any(grepl("H1-hESC", file_name))) {
        sample <- "H1_ESC"  # Standardize to H1_ESC
      } else {
        sample <- paste(parts[(histone_index + 1):(replicate_index - 1)], collapse = "_")
        sample <- standardize_sample_name(sample)  # Apply standardization
      }
      # Calculate widths with explicit conversion
      peaks <- read_delim(file, delim = "\t", col_names = FALSE)
      peaks$X2 <- as.numeric(peaks$X2)
      peaks$X3 <- as.numeric(peaks$X3)
      widths <- peaks$X3 - peaks$X2
      # Append results
      peak_widths <- rbind(peak_widths, 
                           data.frame(Method = method,
                                      Sample = sample,
                                      Histone = histone,
                                      Replicate = replicate,
                                      Width = widths,
                                      DataType = data_type))
    }
  }
  
  return(peak_widths)
}

# Calculate peak widths for both datasets
chip_peak_widths <- calculate_peak_widths(file_paths_chip, "ChIP-seq")
cnr_peak_widths <- calculate_peak_widths(file_paths_cnr, "CUT&RUN")

# Combine both datasets
combined_peak_widths <- rbind(chip_peak_widths, cnr_peak_widths)

combined_peak_widths_filtered <- combined_peak_widths %>%
  filter(
    !Sample %in% c("H1", "H1_Definitive_Endoderm") &
      (!Sample == "K562" | (Sample == "K562" & Replicate %in% c("R1", "R2")))
  ) %>%
  mutate(Sample_Replicate = paste(Sample, Replicate, sep = "_"))

# Apply the same filter to keep only samples that have both data types
combined_peak_widths_filtered <- combined_peak_widths_filtered %>%
  mutate(Sample_Method_Histone_Rep = paste(Sample, Method, Histone, Replicate, sep = "_")) %>%
  filter(Sample_Method_Histone_Rep %in% samples_in_both$Sample_Method_Histone_Rep)

# Create a combined method-datatype variable for faceting
combined_peak_widths_filtered <- combined_peak_widths_filtered %>%
  mutate(Method_DataType = paste(Method, DataType, sep = "_"))

width_violin_plot <- ggplot(combined_peak_widths_filtered, aes(x = Sample_Replicate, y = Width, fill = DataType)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE, scale = "width") +
  facet_grid(Histone ~ Method, scales = "free_x") +
  scale_fill_manual(values = c("ChIP-seq" = "#4C72B0", "CUT&RUN" = "#DD8452")) +
  scale_y_log10(breaks = c(100, 1000, 10000, 100000, 1000000), labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Peak Width Distributions: ChIP-seq vs CUT&RUN",
    x = "Sample_Replicate",
    y = "Peak Width (bp, log10 scale)",
    fill = "Data Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 90),
    panel.spacing = unit(2, "lines")
  ) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.9)) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))

print(width_violin_plot)

# Calculate median peak widths for boxplot comparison
median_widths <- combined_peak_widths_filtered %>%
  group_by(Method, Histone, Sample, Replicate, DataType) %>%
  summarise(MedianWidth = median(Width), .groups = 'drop') %>%
  mutate(Sample_Replicate = paste(Sample, Replicate, sep = "_"))

# Plot median peak widths
median_width_plot <- ggplot(median_widths, aes(x = Sample_Replicate, y = MedianWidth, fill = DataType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.7) +
  facet_grid(Histone ~ Method, scales = "free_x") +
  scale_fill_manual(values = c("ChIP-seq" = "#4C72B0", "CUT&RUN" = "#DD8452")) +
  scale_y_log10(breaks = c(100, 1000, 10000, 100000), labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Median Peak Widths: ChIP-seq vs CUT&RUN",
    x = "Sample_Replicate",
    y = "Median Peak Width (bp, log10 scale)",
    fill = "Data Type"
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
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))

print(median_width_plot)


# Plot ratio of peak counts (CUT&RUN / ChIP-seq)
ratio_data <- combined_peak_counts_merged %>%
  select(Method, Sample_Replicate, Sample_Only, Histone_Only, DataType, Peaks) %>%
  pivot_wider(names_from = DataType, values_from = Peaks) %>%
  mutate(Ratio = `CUT&RUN` / `ChIP-seq`) %>%
  filter(!is.na(Ratio))

ratio_plot <- ggplot(ratio_data, aes(x = Sample_Replicate, y = Ratio, fill = Sample_Only)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.7) +
  facet_grid(Histone_Only ~ Method, scales = "free_x") +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "orange")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Ratio of Peak Counts (CUT&RUN / ChIP-seq)",
    x = "Sample_Replicate",
    y = "Ratio",
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
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))

print(ratio_plot)
dev.off()
