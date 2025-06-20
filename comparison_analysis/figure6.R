library(tidyverse)
library(readxl)
library(ggpubr)

setwd("~/Desktop/R_scripts/")

file_path <- "all_performance_metrics.xlsx"
data <- read_excel(file_path)
data <- data %>% rename(Score = F1)

data_f1 <- data %>% select(Method, Score, HistoneMark) %>% mutate(Metric = "F1")
data_precision <- data %>% select(Method, Precision, HistoneMark) %>% rename(Score = Precision) %>% mutate(Metric = "Precision")
data_recall <- data %>% select(Method, Recall, HistoneMark) %>% rename(Score = Recall) %>% mutate(Metric = "Recall")

plot_metric_by_histone <- function(data, metric_name) {
  ggboxplot(data, x = "Method", y = "Score", fill = "Method", palette = "jco", facet.by = "HistoneMark") +
    stat_compare_means(method = "anova", label = "p.format", label.y = max(data$Score, na.rm = TRUE) * 1.1) +
    stat_compare_means(comparisons = list(c("GOPEAKS", "LANCEOTRON"), c("GOPEAKS", "MACS2"), c("GOPEAKS", "SEACR"),
                                          c("LANCEOTRON", "MACS2"), c("LANCEOTRON", "SEACR"), c("MACS2", "SEACR")),
                       method = "t.test", label = "p.signif", na.rm = TRUE, label.y.npc = "top", label.size = 3) +
    labs(title = paste("Average", metric_name, "Score by Method and Histone Mark"), y = paste(metric_name, "Score"), x = "Method") +
    theme_minimal()
}

plot_f1_histone <- plot_metric_by_histone(data_f1, "F1")
plot_precision_histone <- plot_metric_by_histone(data_precision, "Precision")
plot_recall_histone <- plot_metric_by_histone(data_recall, "Recall")

ggarrange(plot_f1_histone, plot_precision_histone, plot_recall_histone, ncol = 1, nrow = 3)
ggsave("performance_metrics_comparison_by_histone.png", width = 12, height = 18)

