# Load necessary libraries
library(GenomicRanges)
library(readr)
library(ggplot2)

# Define the paths to the BED files
file_paths <- list(
  LANCEOTRON = list(
    "Brain_H3K27ac_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/LANCEOTRON/Brain_H3K27ac_R1_L-tron_noheader.bed",
    "Brain_H3K27ac_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/LANCEOTRON/Brain_H3K27ac_R2_L-tron_noheader.bed",
    "Brain_H3K27me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/LANCEOTRON/Brain_H3K27me3_R1_L-tron_noheader.bed",
    "Brain_H3K27me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/LANCEOTRON/Brain_H3K27me3_R2_L-tron_noheader.bed",
    "Brain_H3K4me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/LANCEOTRON/Brain_H3K4me3_R1_L-tron_noheader.bed",
    "Brain_H3K4me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/LANCEOTRON/Brain_H3K4me3_R2_L-tron_noheader.bed"
  ),
  MACS2 = list(
    "Brain_H3K27ac_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/MACS2/Brain_H3K27ac_R1_peaks.narrowPeak",
    "Brain_H3K27ac_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/MACS2/Brain_H3K27ac_R2_peaks.narrowPeak",
    "Brain_H3K27me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/MACS2/Brain_H3K27me3_R1_peaks.broadPeak",
    "Brain_H3K27me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/MACS2/Brain_H3K27me3_R2_peaks.broadPeak",
    "Brain_H3K4me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/MACS2/Brain_H3K4me3_R1_peaks.narrowPeak",
    "Brain_H3K4me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/MACS2/Brain_H3K4me3_R2_peaks.narrowPeak"
  ),
  GOPEAKS = list(
    "Brain_H3K27ac_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/GOPEAKS/Brain_H3K27ac_R1_peaks_peaks.bed",
    "Brain_H3K27ac_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/GOPEAKS/Brain_H3K27ac_R2_peaks_peaks.bed",
    "Brain_H3K27me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/GOPEAKS/Brain_H3K27me3_R1_peaks_peaks.bed",
    "Brain_H3K27me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/GOPEAKS/Brain_H3K27me3_R2_peaks_peaks.bed",
    "Brain_H3K4me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/GOPEAKS/Brain_H3K4me3_R1_peaks_peaks.bed",
    "Brain_H3K4me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/GOPEAKS/Brain_H3K4me3_R2_peaks_peaks.bed"
  ),
  SEACR = list(
    "Brain_H3K27ac_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/SEACR/Brain_H3K27ac_R1_stringent.stringent.bed",
    "Brain_H3K27ac_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/SEACR/Brain_H3K27ac_R2_stringent.stringent.bed",
    "Brain_H3K27me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/SEACR/Brain_H3K27me3_R1_stringent.stringent.bed",
    "Brain_H3K27me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/SEACR/Brain_H3K27me3_R2_stringent.stringent.bed",
    "Brain_H3K4me3_R1" = "~/Downloads/Transfers/peak_calling_mouse_result/SEACR/Brain_H3K4me3_R1_stringent.stringent.bed",
    "Brain_H3K4me3_R2" = "~/Downloads/Transfers/peak_calling_mouse_result/SEACR/Brain_H3K4me3_R2_stringent.stringent.bed"
  )
)

# Function to read a BED file and count the number of peaks (lines)
count_peaks <- function(file) {
  peaks <- read_delim(file, delim = "\t", col_names = FALSE)
  return(nrow(peaks))
}

# Store peak counts
peak_counts <- data.frame(Method = character(), Replicate = character(), Peaks = integer())

# Loop through files and count peaks
for (method in names(file_paths)) {
  for (replicate in names(file_paths[[method]])) {
    file <- file_paths[[method]][[replicate]]
    num_peaks <- count_peaks(file)
    peak_counts <- rbind(peak_counts, data.frame(Method = method, Replicate = replicate, Peaks = num_peaks))
  }
}


###number of peaks
# Plot the results using ggplot2
# Enhanced bar plot with custom theme and labels
#ggplot(peak_counts, aes(x = Replicate, y = Peaks, fill = Method)) +
#  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.7) + # Adds black border to bars
#  scale_fill_brewer(palette = "Set1") + # Use a colorblind-friendly palette
#  theme_minimal(base_size = 14) + # Set base font size for better readability
#  labs(
#    title = "Peak Counts with Different Methods and Samples",
#    x = "Samples",
#    y = "Number of Peaks",
#    fill = ""
#  ) +
#  theme(
#    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Center and bold the title
#    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), # Rotate x-axis labels for better fit
#    legend.position = "top", # Move the legend to the top
#    panel.grid.major = element_line(color = "gray90") # Lighten the grid lines for a cleaner look
#  ) +
#  geom_text(aes(label = Peaks), position = position_dodge(width = 0.9), vjust = -0.5, size = 4) # Adds peak count labels on top of bars


### Bar plot, R1 and R2 version

# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Assuming peak_counts is your existing data frame
# Modify the data structure
peak_counts_modified <- peak_counts %>%
  separate(Replicate, c("Histone", "Mark", "Replicate"), sep = "_") %>%
  mutate(
    Histone_Mark = paste(Histone, Mark, sep = "_"),
    Method_Replicate = paste(Method, Replicate, sep = "_")
  )

# Create the plot
ggplot(peak_counts_modified, aes(x = Method, y = Peaks, fill = Replicate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.7) +
  facet_wrap(~Histone_Mark, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("R1" = "#66c2a5", "R2" = "#fc8d62")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Peak Counts with Different Methods and Samples",
    x = "Method",
    y = "Number of Peaks",
    fill = "Replicate"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  ) +
  geom_text(aes(label = Peaks), position = position_dodge(width = 0.9), 
            vjust = -0.5, hjust = 0.5, size = 3, angle = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))





###peak width

# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

# Function to calculate peak widths from BED files
calculate_peak_widths <- function(file) {
  peaks <- read_delim(file, delim = "\t", col_names = FALSE)
  widths <- peaks$X3 - peaks$X2
  return(widths)  # Return all widths
}

# Store all peak widths
peak_widths <- data.frame(Method = character(), Replicate = character(), Width = numeric())

# Loop through files to calculate peak widths
for (method in names(file_paths)) {
  for (replicate in names(file_paths[[method]])) {
    file <- file_paths[[method]][[replicate]]
    widths <- calculate_peak_widths(file)
    peak_widths <- rbind(peak_widths, data.frame(Method = method, Replicate = replicate, Width = widths))
  }
}

# Modify the data structure
peak_widths_modified <- peak_widths %>%
  separate(Replicate, c("Histone", "Mark", "Replicate"), sep = "_") %>%
  mutate(
    Histone_Mark = paste(Histone, Mark, sep = "_"),
    Method_Replicate = paste(Method, Replicate, sep = "_")
  )

# Create the plot
ggplot(peak_widths_modified, aes(x = Method, y = Width, fill = Replicate)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE, scale = "width") +
  facet_wrap(~Histone_Mark, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("R1" = "#66c2a5", "R2" = "#fc8d62")) +
  scale_y_log10() +  # Use log scale for better visualization of differences
  theme_minimal(base_size = 14) +
  labs(
    title = "Log-ScaledPeak Width Distributions with Different Methods and Samples",
    x = "Method",
    y = "Peak Width (log10 bp)",
    fill = "Replicate"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  ) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.9)) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))

# Save the plot
#ggsave("peak_width_distribution_mouse.pdf", width = 13, height = 8)




##consensu peaks
