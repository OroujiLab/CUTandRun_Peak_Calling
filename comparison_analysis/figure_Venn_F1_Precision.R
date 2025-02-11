library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(ggVennDiagram)
library(viridis)
library(gridExtra)

file_paths <- list(
  LANCEOTRON = Sys.glob("~/Downloads/Transfers/results_2/LANCEOTRON/*H3K*_R*.bed"),
  MACS2 = c(
    Sys.glob("~/Downloads/Transfers/results_2/MACS2/*H3K27me3*_R*.broadPeak"),
    Sys.glob("~/Downloads/Transfers/results_2/MACS2/*H3K*_R*.narrowPeak")
  ),
  GOPEAKS = Sys.glob("~/Downloads/Transfers/results_2/GOPEAKS/*H3K*_R*_peaks.bed"),
  SEACR = Sys.glob("~/Downloads/Transfers/results_2/SEACR/*H3K*_R*_stringent*.bed")
)

get_sample_info <- function(filename) {
  tryCatch({
    parts <- strsplit(basename(filename), "_")[[1]]
    id <- parts[1]  
    histone_idx <- grep("H3K", parts)[1]
    replicate_idx <- grep("^R\\d+", parts)[1]
    
    if (is.na(histone_idx) || is.na(replicate_idx)) {
      stop("Invalid filename format")
    }
    
    histone <- parts[histone_idx]  
    sample_name <- paste(parts[(histone_idx + 1):(replicate_idx - 1)], collapse = "_")  
    replicate <- parts[replicate_idx]  
    
    return(list(
      ID = id,
      HistoneMark = histone,
      Sample = sample_name,
      Replicate = replicate
    ))
  }, error = function(e) {
    warning(sprintf("Error parsing filename %s: %s", filename, e$message))
    return(list(ID = NA_character_, 
                HistoneMark = NA_character_,
                Sample = NA_character_, 
                Replicate = NA_character_))
  })
}

get_sample_name <- function(filename) {
  tryCatch({
    parts <- strsplit(basename(filename), "_")[[1]]
    histone_idx <- grep("H3K", parts)[1]
    replicate_idx <- grep("^R\\d+", parts)[1]
    
    if (is.na(histone_idx) || is.na(replicate_idx)) {
      return(NA_character_)
    }
    
    sample_name <- paste(parts[(histone_idx + 1):(replicate_idx - 1)], collapse = "_")
    return(sample_name)
  }, error = function(e) {
    warning(sprintf("Error getting sample name from %s: %s", filename, e$message))
    return(NA_character_)
  })
}

read_peaks <- function(file) {
  peaks <- read_delim(file, delim = "\t", col_names = FALSE, show_col_types = FALSE)
  peaks_subset <- data.frame(
    chr = peaks$X1,
    start = peaks$X2,
    end = peaks$X3
  )
  return(peaks_subset)
}

# Function to write BED file
write_bed <- function(gr, filename) {
  df <- data.frame(
    chr = seqnames(gr),
    start = start(gr),
    end = end(gr)
  )
  write.table(df, filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Function to create Venn diagram and save overlaps
create_venn_diagram <- function(method_peaks, histone_mark, sample_name, output_dir) {
  dir.create(file.path(output_dir, sample_name, histone_mark), recursive = TRUE, showWarnings = FALSE)
  
  # Convert all peaks to GRanges objects
  peak_granges <- lapply(method_peaks, function(peaks) {
    GRanges(seqnames = peaks$chr,
            ranges = IRanges(start = peaks$start, end = peaks$end))
  })
  
  # Get all possible 3-method combinations
  methods <- names(peak_granges)
  combinations <- combn(methods, 3)
  combination_consensuses <- GRangesList()
  
  # Get consensus for each combination
  for(i in 1:ncol(combinations)) {
    combo <- combinations[,i]
    combo_name <- paste(combo, collapse="_")
    
    # Get intersection of all three methods
    peaks1 <- peak_granges[[combo[1]]]
    peaks2 <- peak_granges[[combo[2]]]
    peaks3 <- peak_granges[[combo[3]]]
    
    # Find regions that overlap in all three methods
    overlap12 <- subsetByOverlaps(peaks1, peaks2)
    combo_consensus <- subsetByOverlaps(overlap12, peaks3)
    
    if(length(combo_consensus) > 0) {
      combination_consensuses[[combo_name]] <- combo_consensus
      
      # Save individual combination consensus peaks
      combo_file <- file.path(output_dir, sample_name, histone_mark, 
                              paste0("consensus_", combo_name, ".bed"))
      write_bed(combo_consensus, combo_file)
    }
  }
  
  # Create final consensus as union of all combination consensuses
  if(length(combination_consensuses) > 0) {
    all_consensus_peaks <- unlist(combination_consensuses)
    consensus_regions <- reduce(all_consensus_peaks)
  } else {
    consensus_regions <- GRanges()
  }
  # Find overlaps for each method
  venn_lists <- list()
  overlap_counts <- list()
  
  # Calculate metrics for each method
  metrics <- data.frame()
  
  for(method in names(peak_granges)) {
    if(!is.null(peak_granges[[method]])) {
      # For Venn diagram
      overlaps <- findOverlaps(peak_granges[[method]], consensus_regions)
      venn_lists[[method]] <- unique(subjectHits(overlaps))
      
      # Calculate metrics
      true_positives <- length(unique(queryHits(overlaps)))
      false_positives <- length(peak_granges[[method]]) - true_positives
      false_negatives <- length(consensus_regions) - length(unique(subjectHits(overlaps)))
      
      precision <- true_positives / (true_positives + false_positives)
      recall <- true_positives / (true_positives + false_negatives)
      f1 <- 2 * (precision * recall) / (precision + recall)
      
      metrics <- rbind(metrics, data.frame(
        Method = method,
        Sample = sample_name,
        HistoneMark = histone_mark,
        Precision = round(precision, 3),
        Recall = round(recall, 3),
        F1 = round(f1, 3),
        TP = true_positives,
        FP = false_positives,
        FN = false_negatives
      ))
      
      # Save individual method peaks
      out_file <- file.path(output_dir, sample_name, histone_mark, 
                            paste0(method, "_peaks.bed"))
      write_bed(peak_granges[[method]], out_file)
    }
  }
  
  # Find and save overlapping peaks between pairs
  for(method1 in names(peak_granges)) {
    for(method2 in names(peak_granges)) {
      if(method1 < method2) {  # Avoid duplicate comparisons
        overlap <- GenomicRanges::reduce(subsetByOverlaps(peak_granges[[method1]], 
                                                          peak_granges[[method2]]))
        
        # Save overlapping peaks
        out_file <- file.path(output_dir, sample_name, histone_mark,
                              paste0(method1, "_", method2, "_overlap.bed"))
        write_bed(overlap, out_file)
        
        # Store overlap count
        overlap_counts[[paste(method1, method2, sep="_")]] <- length(overlap)
      }
    }
  }
  
  cat("\nConsensus Statistics for each 3-method combination:\n")
  combo_stats <- data.frame()
  for(combo_name in names(combination_consensuses)) {
    n_peaks <- length(combination_consensuses[[combo_name]])
    cat(sprintf("%s: %d peaks\n", combo_name, n_peaks))
    combo_stats <- rbind(combo_stats, data.frame(
      Sample = sample_name,
      HistoneMark = histone_mark,
      Combination = combo_name,
      NumPeaks = n_peaks
    ))
  }
  
  combo_stats_file <- file.path(output_dir, sample_name, histone_mark, "combination_stats.txt")
  write.table(combo_stats, combo_stats_file, sep="\t", quote=FALSE, row.names=FALSE)
  
  consensus_file <- file.path(output_dir, sample_name, histone_mark, "consensus_peaks.bed")
  write_bed(consensus_regions, consensus_file)
  
  # Print overlap statistics
  cat("\nOverlap Statistics for", histone_mark, "-", sample_name, ":\n")
  for(comparison in names(overlap_counts)) {
    cat(sprintf("%s: %d peaks\n", comparison, overlap_counts[[comparison]]))
  }
  
  cat("\nPerformance Metrics (using consensus peaks as ground truth):\n")
  print(metrics)
  
  metrics_file <- file.path(output_dir, sample_name, histone_mark, "performance_metrics.txt")
  write.table(metrics, metrics_file, sep="\t", quote=FALSE, row.names=FALSE)
  
  venn_plot <- ggVennDiagram(venn_lists) +
    scale_fill_viridis() +
    labs(title = paste(histone_mark, "-", sample_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Create performance metrics plot
  metrics_long <- tidyr::pivot_longer(metrics, 
                                      cols = c("Precision", "Recall", "F1"),
                                      names_to = "Metric", 
                                      values_to = "Value")
  
  metrics_plot <- ggplot(metrics_long, aes(x = Method, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis_d() +
    ylim(0, 1) +
    labs(title = paste("Performance Metrics -", histone_mark, "-", sample_name),
         y = "Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf_file <- file.path(output_dir, sample_name, histone_mark, "venn_diagram.pdf")
  pdf(pdf_file, width = 10, height = 8)
  print(venn_plot)
  dev.off()
  
  perf_file <- file.path(output_dir, sample_name, histone_mark, "performance_metrics.pdf")
  pdf(perf_file, width = 8, height = 6)
  print(metrics_plot)
  dev.off()
  
  return(list(
    venn = venn_plot, 
    metrics = metrics_plot,
    performance_metrics = metrics,
    combination_stats = combo_stats
  ))
}

analyze_histone_mark_sample <- function(file_paths, histone_mark, sample_name, output_dir) {
  method_peaks <- list()
  
  for (method in c("LANCEOTRON", "MACS2", "GOPEAKS", "SEACR")) {
    files <- file_paths[[method]]
    mark_sample_files <- files[grep(paste0(histone_mark, ".*", sample_name), files)]
    
    if (length(mark_sample_files) > 0) {
      peaks_data <- data.frame()
      for (file in mark_sample_files) {
        peaks <- read_peaks(file)
        peaks_data <- rbind(peaks_data, peaks)
      }
      if (nrow(peaks_data) > 0) {
        method_peaks[[method]] <- unique(peaks_data)
        cat(sprintf("%s - %s peaks: %d\n", method, histone_mark, nrow(peaks_data)))
      }
    }
  }
  
  if (length(method_peaks) < 2) {
    cat(sprintf("\nNot enough methods with peaks for %s - %s\n", histone_mark, sample_name))
    return(NULL)
  }
  
  return(create_venn_diagram(method_peaks, histone_mark, sample_name, output_dir))
}

# Main function to generate all analyses
generate_peak_analyses <- function(file_paths, output_dir = "peak_analyses") {
  histone_marks <- c("H3K27ac", "H3K27me3", "H3K4me3")
  samples <- unique(sapply(unlist(file_paths), get_sample_name))
  samples <- samples[!is.na(samples)]
  
  cat("Processing the following samples:", paste(samples, collapse=", "), "\n")
  cat("For histone marks:", paste(histone_marks, collapse=", "), "\n\n")
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # Store all plots and metrics
  all_venn_plots <- list()
  all_metrics_plots <- list()
  all_performance_metrics <- data.frame()
  all_combination_stats <- data.frame()
  
  for (sample in samples) {
    cat(sprintf("\nProcessing sample: %s\n", sample))
    for (mark in histone_marks) {
      cat(sprintf("\nProcessing mark: %s\n", mark))
      results <- analyze_histone_mark_sample(file_paths, mark, sample, output_dir)
      if (!is.null(results)) {
        all_venn_plots[[paste(sample, mark, sep="_")]] <- results$venn
        all_metrics_plots[[paste(sample, mark, sep="_")]] <- results$metrics
        all_performance_metrics <- rbind(all_performance_metrics, results$performance_metrics)
        all_combination_stats <- rbind(all_combination_stats, results$combination_stats)
      }
    }
  }
  
  write.csv(all_performance_metrics, 
            file.path(output_dir, "all_performance_metrics.csv"), 
            row.names = FALSE)
  write.csv(all_combination_stats, 
            file.path(output_dir, "all_combination_stats.csv"), 
            row.names = FALSE)
  
  # Save all Venn diagrams
  if (length(all_venn_plots) > 0) {
    pdf(file.path(output_dir, "all_venn_diagrams.pdf"), width = 10, height = 8)
    for (plot in all_venn_plots) {
      print(plot)
    }
    dev.off()
  }
  
  # Save all metrics plots
  if (length(all_metrics_plots) > 0) {
    pdf(file.path(output_dir, "all_performance_metrics.pdf"), width = 10, height = 6)
    for (plot in all_metrics_plots) {
      print(plot)
    }
    dev.off()
  }
  
  cat("\nAnalyses completed. Results saved in:", output_dir, "\n")
  #cat("\nGenerated files:\n")
  #cat("1. all_performance_metrics.csv - Combined performance metrics for all samples\n")
  #cat("2. all_combination_stats.csv - Statistics for all 3-method combinations\n")
  #cat("3. all_venn_diagrams.pdf - Combined Venn diagrams\n")
  #cat("4. all_performance_metrics.pdf - Combined performance plots\n")
}

generate_peak_analyses(file_paths)
