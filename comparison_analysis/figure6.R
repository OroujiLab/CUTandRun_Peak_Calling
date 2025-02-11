# Benchmarking Peak Calling Method for CUT&RUN
# Author: Elias Orouji
# Date: November 9, 2024
# Manuscript Title: "Benchmarking Peak Calling Methods for CUT&RUN"
# Figure Number: [Figure 6]

# Description: This script generates ternary plots for precision, recall and F1 score for all the four peak calling methods based on the histone mark tested.

# Notes:
# - This script is intended for use with data from mouse brain samples with three different histone marks.
# - The peak calling methods compared include MACS2, SEACR, GoPeaks, and LanceOtron.
# - Ensure the appropriate input files are loaded for each figure, and customize the parameters as needed.

# Clean up environment
rm(list = ls())

# Install required packages if they are not already installed
if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

# Load necessary libraries
library(ggtern)
library(gridExtra)

# Create data frames for each plot
data_precision <- data.frame(
  H3K27ac = c(0.316741, 0.314521, 0.282387, 0.170911, 0.308261, 0.198829, 0.269176, 0.194383, 0.288525, 0.314426, 0.267023, 0.183279, 0.301659, 0.292911, 0.252298, 0.178029, 0.269618, 0.220763, 0.290771, 0.261804),
  H3K4me3 = c(0.286303, 0.184429, 0.282965, 0.282965, 0.298705, 0.173790, 0.275330, 0.272852, 0.313443, 0.244059, 0.272331, 0.272059, 0.272398, 0.244893, 0.276014, 0.276014, 0.286107, 0.290139, 0.282160, 0.185751),
  H3K27me3 = c(0.082435, 0.228462, 0.304722, 0.263158, 0.194205, 0.273988, 0.283046, 0.257434, 0.083607, 0.267023, 0.221896, 0.272331, 0.133032, 0.247957, 0.254852, 0.269942, 0.223511, 0.116941, 0.302149, 0.270286),
  Color = rep(c("#1f78b4", "#ffcc00", "#4d4d4d", "#e31a1c"), 5),
  Label = rep(c("GoPeaks", "LanceOtron", "MACS2", "SEACR"), 5)
)

data_recall <- data.frame(
  H3K27ac = c(0.103495, 0.335349, 0.325269, 0.235887, 0.040362, 0.359834, 0.313775, 0.165849, 0.286343, 0.248980, 0.316444, 0.206664, 0.208869, 0.154976, 0.203695, 0.166731, 0.245687, 0.245617, 0.304563, 0.309234),
  H3K4me3 = c(0.172098, 0.366166, 0.365800, 0.095936, 0.192120, 0.370841, 0.252782, 0.322300, 0.284690, 0.268284, 0.286206, 0.302429, 0.303508, 0.276620, 0.268299, 0.297748, 0.269124, 0.284394, 0.214789, 0.274623),
  H3K27me3 = c(0.242965, 0.342141, 0.195951, 0.218943, 0.301851, 0.269894, 0.199022, 0.207762, 0.186146, 0.188920, 0.170161, 0.235592, 0.292643, 0.234901, 0.225581, 0.220998, 0.247408, 0.229989, 0.480648, 0.416143),
  Color = rep(c("#1f78b4", "#ffcc00", "#4d4d4d", "#e31a1c"), 5),
  Label = rep(c("GoPeaks", "LanceOtron", "MACS2", "SEACR"), 5)
)

data_f1 <- data.frame(
  H3K27ac = c(0.178476696, 0.156119742, 0.35392194, 0.311481622, 0.078138718, 0.339332748, 0.247585601, 0.334942932, 0.350134048, 0.218230563, 0.038605898, 0.393029491, 0.115102356, 0.235998455, 0.3198146, 0.329084589, 0.291826309, 0.261813538, 0.179438059, 0.266922095),
  H3K4me3 = c(0.230685921, 0.271841155, 0.34765343, 0.149819495, 0.229420505, 0.24576523, 0.274294205, 0.250520059, 0.270932607, 0.273995916, 0.184819605, 0.270251872, 0.26252505, 0.261857047, 0.308283233, 0.167334669, 0.260639777, 0.25674548, 0.262030598, 0.220584145),
  H3K27me3 = c(0.24526352, 0.252152945, 0.239063038, 0.263520496, 0.338491296, 0.294003868, 0.33926499, 0.028239845, 0.201673272, 0.399383531, 0, 0.398943197, 0.242726081, 0.255832241, 0.245609436, 0.255832241, 0.205089252, 0.307633878, 0.126471705, 0.360805165),
  Color = rep(c("#1f78b4", "#ffcc00", "#4d4d4d", "#e31a1c"), 5),
  Label = rep(c("GoPeaks", "LanceOtron", "MACS2", "SEACR"), 5)
)

# Function to create ternary plots
create_ternary_plot <- function(data, title) {
  ggtern(data = data, aes(x = H3K27ac, y = H3K4me3, z = H3K27me3)) +
    geom_point(aes(color = Color), size = 5, alpha = 0.8) +
    scale_color_identity(name = "Method", breaks = unique(data$Color), labels = unique(data$Label), guide = "legend") +
    labs(
      title = title,
      x = "H3K27ac",
      y = "H3K4me3",
      z = "H3K27me3"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    theme_showgrid()
}

# Generate each plot individually
plot1 <- create_ternary_plot(data_precision, "Peak Caller's Precision for Each Histone Mark")
print(plot1)

plot2 <- create_ternary_plot(data_recall, "Peak Caller's Recall for Each Histone Mark")
print(plot2)

plot3 <- create_ternary_plot(data_f1, "Peak Caller's F1 Score for Each Histone Mark")
print(plot3)

