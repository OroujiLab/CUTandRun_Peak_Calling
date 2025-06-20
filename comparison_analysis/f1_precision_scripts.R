### versioni ke avalin valuesh ba 0.63 shoroo mshe

library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(viridis)


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
    id <- parts[1]  # Gets "4DNFI258RO3L"
    histone_idx <- grep("H3K", parts)[1]
    replicate_idx <- grep("^R\\d+", parts)[1]
    
    if (is.na(histone_idx) || is.na(replicate_idx)) {
      stop("Invalid filename format")
    }
    
    histone <- parts[histone_idx]  # Gets "H3K27ac"
    sample_name <- paste(parts[(histone_idx + 1):(replicate_idx - 1)], collapse = "_")  # Gets "H1_ESC"
    replicate <- parts[replicate_idx]  # Gets "R1"
    
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

write_bed <- function(gr, filename) {
  df <- data.frame(
    chr = seqnames(gr),
    start = start(gr),
    end = end(gr)
  )
  write.table(df, filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


calculate_performance_metrics <- function(peak_granges, all_regions) {
  # Create consensus peaks (regions found by at least 3 methods)
  peak_coverage <- coverage(unlist(GRangesList(peak_granges)))
  consensus_peaks <- GenomicRanges::reduce(all_regions[peak_coverage >= 3])
  
  metrics <- data.frame()
  for(method in names(peak_granges)) {
    # Find overlaps with consensus
    overlaps <- findOverlaps(peak_granges[[method]], consensus_peaks)
    
    true_positives <- length(unique(queryHits(overlaps)))
    false_positives <- length(peak_granges[[method]]) - true_positives
    false_negatives <- length(consensus_peaks) - length(unique(subjectHits(overlaps)))
    
    precision <- true_positives / (true_positives + false_positives)
    recall <- true_positives / (true_positives + false_negatives)
    f1 <- 2 * (precision * recall) / (precision + recall)
    
    metrics <- rbind(metrics, data.frame(
      Method = method,
      Precision = round(precision, 3),
      Recall = round(recall, 3),
      F1 = round(f1, 3)
    ))
  }
  return(list(metrics = metrics, consensus_peaks = consensus_peaks))
}
create_venn_diagram <- function(method_peaks, histone_mark, sample_name, output_dir) {
  dir.create(file.path(output_dir, sample_name, histone_mark), recursive = TRUE, showWarnings = FALSE)
  
  peak_granges <- lapply(method_peaks, function(peaks) {
    GRanges(seqnames = peaks$chr,
            ranges = IRanges(start = peaks$start, end = peaks$end))
  })
  
  all_peaks <- unlist(GRangesList(peak_granges))
  all_regions <- GenomicRanges::reduce(all_peaks)
  
  cvg <- coverage(all_peaks)
  consensus_regions <- GRanges()
  
  for(seqname in names(cvg)) {
    runs <- slice(cvg[[seqname]], lower=3)
    if(length(runs@ranges) > 0) {
      consensus_regions <- c(consensus_regions,
                             GRanges(seqname,
                                     runs@ranges))
    }
  }
  
  venn_lists <- list()
  overlap_counts <- list()
  
  metrics <- data.frame()
  
  for(method in names(peak_granges)) {
    if(!is.null(peak_granges[[method]])) {
      overlaps <- findOverlaps(peak_granges[[method]], all_regions)
      venn_lists[[method]] <- unique(subjectHits(overlaps))
      
      if(length(consensus_regions) > 0) {
        consensus_overlaps <- findOverlaps(peak_granges[[method]], consensus_regions)
        
        true_positives <- length(unique(queryHits(consensus_overlaps)))
        false_positives <- length(peak_granges[[method]]) - true_positives
        false_negatives <- length(consensus_regions) - length(unique(subjectHits(consensus_overlaps)))
        
        precision <- true_positives / (true_positives + false_positives)
        recall <- true_positives / (true_positives + false_negatives)
        f1 <- 2 * (precision * recall) / (precision + recall)
        
        metrics <- rbind(metrics, data.frame(
          Method = method,
          Precision = round(precision, 3),
          Recall = round(recall, 3),
          F1 = round(f1, 3),
          TP = true_positives,
          FP = false_positives,
          FN = false_negatives
        ))
      }
      
      out_file <- file.path(output_dir, sample_name, histone_mark, 
                            paste0(method, "_peaks.bed"))
      write_bed(peak_granges[[method]], out_file)
    }
  }
  
  for(method1 in names(peak_granges)) {
    for(method2 in names(peak_granges)) {
      if(method1 < method2) {  # Avoid duplicate comparisons
        overlap <- GenomicRanges::reduce(subsetByOverlaps(peak_granges[[method1]], 
                                                          peak_granges[[method2]]))
        
        out_file <- file.path(output_dir, sample_name, histone_mark,
                              paste0(method1, "_", method2, "_overlap.bed"))
        write_bed(overlap, out_file)
        
        overlap_counts[[paste(method1, method2, sep="_")]] <- length(overlap)
      }
    }
  }
  
  if(length(consensus_regions) > 0) {
    consensus_file <- file.path(output_dir, sample_name, histone_mark, "consensus_peaks.bed")
    write_bed(consensus_regions, consensus_file)
  }
  
  cat("\nOverlap Statistics for", histone_mark, "-", sample_name, ":\n")
  for(comparison in names(overlap_counts)) {
    cat(sprintf("%s: %d peaks\n", comparison, overlap_counts[[comparison]]))
  }
  
  if(nrow(metrics) > 0) {
    cat("\nPerformance Metrics (using consensus peaks as ground truth):\n")
    print(metrics)
    
    metrics_file <- file.path(output_dir, sample_name, histone_mark, "performance_metrics.txt")
    write.table(metrics, metrics_file, sep="\t", quote=FALSE, row.names=FALSE)
  }
  
  venn_plot <- ggVennDiagram(venn_lists) +
    scale_fill_viridis() +
    labs(title = paste(histone_mark, "-", sample_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(nrow(metrics) > 0) {
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
    
    perf_file <- file.path(output_dir, sample_name, histone_mark, "performance_metrics.pdf")
    pdf(perf_file, width = 8, height = 6)
    print(metrics_plot)
    dev.off()
  } else {
    metrics_plot <- NULL
  }
  
  pdf_file <- file.path(output_dir, sample_name, histone_mark, "venn_diagram.pdf")
  pdf(pdf_file, width = 10, height = 8)
  print(venn_plot)
  dev.off()
  
  return(list(venn = venn_plot, metrics = metrics_plot))
}



generate_peak_analyses <- function(file_paths, output_dir = "peak_analyses") {
  histone_marks <- c("H3K27ac", "H3K27me3", "H3K4me3")
  samples <- unique(sapply(unlist(file_paths), get_sample_name))
  samples <- samples[!is.na(samples)]
  
  cat("Processing the following samples:", paste(samples, collapse=", "), "\n")
  cat("For histone marks:", paste(histone_marks, collapse=", "), "\n\n")
  
  dir.create(output_dir, showWarnings = FALSE)
  all_venn_plots <- list()
  all_metrics_plots <- list()
  
  for (sample in samples) {
    cat(sprintf("\nProcessing sample: %s\n", sample))
    for (mark in histone_marks) {
      cat(sprintf("\nProcessing mark: %s\n", mark))
      plots <- analyze_histone_mark_sample(file_paths, mark, sample, output_dir)
      if (!is.null(plots)) {
        all_venn_plots[[paste(sample, mark, sep="_")]] <- plots$venn
        all_metrics_plots[[paste(sample, mark, sep="_")]] <- plots$metrics
      }
    }
  }
  
  if (length(all_venn_plots) > 0) {
    pdf(file.path(output_dir, "all_venn_diagrams.pdf"), width = 10, height = 8)
    for (plot in all_venn_plots) {
      print(plot)
    }
    dev.off()
  }
  
  if (length(all_metrics_plots) > 0) {
    pdf(file.path(output_dir, "all_performance_metrics.pdf"), width = 10, height = 6)
    for (plot in all_metrics_plots) {
      print(plot)
    }
    dev.off()
  }
  
  all_metrics <- data.frame()
  
  for (sample in samples) {
    for (mark in histone_marks) {
      metrics_file <- file.path(output_dir, sample, mark, "performance_metrics.txt")
      if(file.exists(metrics_file)) {
        metrics <- read.table(metrics_file, header=TRUE, sep="\t")
        metrics$Sample <- sample
        metrics$HistoneMark <- mark
        all_metrics <- rbind(all_metrics, metrics)
      }
    }
  }
  all_metrics <- all_metrics[, c("Sample", "HistoneMark", "Method", 
                                 "Precision", "Recall", "F1", "TP", "FP", "FN")]
  
  write.csv(all_metrics, 
            file.path(output_dir, "all_performance_metrics.csv"), 
            row.names = FALSE)
  
  cat("\nComprehensive metrics saved to: all_performance_metrics.csv\n")}

generate_peak_analyses(file_paths)