library(tidyverse)
library(GenomicRanges)
library(ComplexUpset)
library(dplyr)
library(ggplot2)


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
    id <- parts[1]  # gets first part of filenames
    histone_idx <- grep("H3K", parts)[1]
    replicate_idx <- grep("^R\\d+", parts)[1]
    
    if (is.na(histone_idx) || is.na(replicate_idx)) {
      stop("Invalid filename format")
    }
    
    histone <- parts[histone_idx]  # Gets the mark
    sample_name <- paste(parts[(histone_idx + 1):(replicate_idx - 1)], collapse = "_")  # Gets "H1_ESC"
    replicate <- parts[replicate_idx]  # Gets the replicate
    
    return(list(
      ID = id,
      HistoneMark = histone,
      Sample = sample_name,
      Replicate = replicate
    ))
  }, error = function(e) {
    warning(sprintf("Error filename %s: %s", filename, e$message))
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
    
    # Get parts between histone mark and replicate number which is the sample name
    sample_name <- paste(parts[(histone_idx + 1):(replicate_idx - 1)], collapse = "_")
    return(sample_name)
  }, error = function(e) {
    warning(sprintf("Error getting sample name from %s: %s", filename, e$message))
    return(NA_character_)
  })
}

# reading peak files, only first 3 columns
read_peaks <- function(file) {
  peaks <- read_delim(file, delim = "\t", col_names = FALSE, show_col_types = FALSE)
  peaks_subset <- data.frame(
    chr = peaks$X1,
    start = peaks$X2,
    end = peaks$X3
  )
  return(peaks_subset)
}

df_to_granges <- function(df) {
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end)
  )
}


#test_name_extraction <- function(filename) {
#  info <- get_sample_info(filename)
#  cat("File:", filename, "\n")
#  cat("ID:", info$ID, "\n")
#  cat("Histone Mark:", info$HistoneMark, "\n")
#  cat("Sample:", info$Sample, "\n")
#  cat("Replicate:", info$Replicate, "\n\n")
#  
#  # Also test get_sample_name
#  sample <- get_sample_name(filename)
#  cat("get_sample_name result:", sample, "\n")
#  cat("------------------\n")
#}

analyze_histone_mark_sample <- function(file_paths, histone_mark, sample_name) {
  results <- list()
  method_peaks <- list()
  method_peak_counts <- list()
  for (method in names(file_paths)) {
    files <- file_paths[[method]]
    # Modified pattern to match ID_histone_sample_R format
    mark_sample_files <- files[grep(paste0("H3K", "[^_]+_", sample_name, "_R"), files)]
    
    if (length(mark_sample_files) > 0) {
      peaks_data <- data.frame()
      for (file in mark_sample_files) {
        peaks <- read_peaks(file)
        peaks_data <- rbind(peaks_data, peaks)
      }
      if (nrow(peaks_data) > 0) {
        peaks_data <- unique(peaks_data)
        method_peaks[[method]] <- peaks_data
        method_peak_counts[[method]] <- nrow(peaks_data)
      }
    }
  }
  
  method_peaks <- method_peaks[sapply(method_peaks, nrow) > 0]
  methods <- names(method_peaks)
  
  if (length(methods) < 2) {
    cat(sprintf("\nNot enough methods with peaks for %s - %s\n", histone_mark, sample_name))
    return(NULL)
  }
  
  for (method in methods) {
    current_df <- method_peaks[[method]]
    other_peaks <- do.call(rbind, method_peaks[names(method_peaks) != method])
    
    if (!is.null(other_peaks) && nrow(other_peaks) > 0) {
      current_gr <- df_to_granges(current_df)
      other_gr <- df_to_granges(other_peaks)
      
      overlaps <- findOverlaps(current_gr, other_gr)
      unique_indices <- setdiff(seq_len(nrow(current_df)), queryHits(overlaps))
      unique_peaks <- current_df[unique_indices, ]
    } else {
      unique_peaks <- current_df
    }
    
    results[[paste0(method, "_unique")]] <- nrow(unique_peaks)
  }
  for (i in 1:(length(methods)-1)) {
    for (j in (i+1):length(methods)) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      
      overlaps <- findOverlaps(gr1, gr2)
      overlap_count <- length(unique(queryHits(overlaps)))
      
      results[[paste(method1, method2, "shared", sep = "_")]] <- overlap_count
      cat(sprintf("%s-%s shared peaks: %d\n", method1, method2, overlap_count))
    }
  }
  if (length(methods) >= 3) {
    cat("Finding three-way overlaps...\n")
    combinations <- combn(methods, 3, simplify = FALSE)
    
    for (combo in combinations) {
      method1 <- combo[1]
      method2 <- combo[2]
      method3 <- combo[3]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      gr3 <- df_to_granges(method_peaks[[method3]])
      
      overlaps_12 <- findOverlaps(gr1, gr2)
      shared_gr_12 <- gr1[unique(queryHits(overlaps_12))]
      
      overlaps_123 <- findOverlaps(shared_gr_12, gr3)
      overlap_count <- length(unique(queryHits(overlaps_123)))
      
      key <- paste(method1, method2, method3, "shared", sep = "_")
      results[[key]] <- overlap_count
      cat(sprintf("%s-%s-%s shared peaks: %d\n", method1, method2, method3, overlap_count))
    }
  }
  
  if (length(methods) >= 2) {
    cat("Finding consensus peaks...\n")
    gr_list <- lapply(methods, function(m) df_to_granges(method_peaks[[m]]))
    consensus_gr <- gr_list[[1]]
    for (i in 2:length(gr_list)) {
      overlaps <- findOverlaps(consensus_gr, gr_list[[i]])
      consensus_gr <- consensus_gr[unique(queryHits(overlaps))]
      
      if (length(consensus_gr) == 0) break
    }
    
    results[["consensus"]] <- length(consensus_gr)
    cat(sprintf("Consensus peaks (shared among all methods): %d\n", length(consensus_gr)))
  }
  
  generate_peak_files(method_peaks, results, histone_mark, sample_name)
  csv_file_path <- sprintf("%s_%s_peaks_summary.csv", sample_name, histone_mark)
  write_results_to_csv(results, method_peak_counts, histone_mark, sample_name, csv_file_path)
  
  return(list(results = results, method_peak_counts = method_peak_counts))
}


check_available_files <- function(file_paths) {
  cat("\nChecking available files:\n")
  for (method in names(file_paths)) {
    files <- file_paths[[method]]
    cat(sprintf("\nMethod: %s\n", method))
    cat("Files found:\n")
    for (file in files) {
      cat(sprintf("- %s\n", basename(file)))
    }
  }
  all_files <- unlist(file_paths)
  all_bases <- sapply(all_files, function(f) {
    info <- get_sample_info(f)
    paste(info$HistoneMark, info$Sample, info$Replicate, sep = "_")
  })
  
  duplicates <- all_bases[duplicated(all_bases)]
  if (length(duplicates) > 0) {
    cat("\nWarning: Potential duplicate files found:\n")
    for (dup in unique(duplicates)) {
      cat(sprintf("- %s\n", dup))
    }
  }
}

write_peaks_to_bed <- function(peaks_df, file_path) {
  write.table(peaks_df, file = file_path, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

generate_peak_files <- function(method_peaks, results, histone_mark, sample_name) {
  dir.create("peak_files", showWarnings = FALSE)
  methods <- names(method_peaks)
  granges_to_df <- function(gr) {
    data.frame(
      chr = seqnames(gr),
      start = start(gr),
      end = end(gr)
    )
  }
  
  #Generate unique peak files 
  for (method in methods) {
    current_df <- method_peaks[[method]]
    other_peaks <- do.call(rbind, method_peaks[names(method_peaks) != method])
    
    current_gr <- df_to_granges(current_df)
    if (!is.null(other_peaks) && nrow(other_peaks) > 0) {
      other_gr <- df_to_granges(other_peaks)
      overlaps <- findOverlaps(current_gr, other_gr)
      unique_indices <- setdiff(seq_len(length(current_gr)), queryHits(overlaps))
      unique_peaks <- granges_to_df(current_gr[unique_indices])
    } else {
      unique_peaks <- current_df
    }
    
    file_path <- sprintf("peak_files/%s_%s_%s_unique_peaks_all.bed", 
                         sample_name, histone_mark, method)
    write_peaks_to_bed(unique_peaks, file_path)
    cat(sprintf("Written %d unique peaks for %s to %s\n", 
                nrow(unique_peaks), method, file_path))
  }
  
  #Generate pairwise shared peak file
  for (i in 1:(length(methods)-1)) {
    for (j in (i+1):length(methods)) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      
      overlaps <- findOverlaps(gr1, gr2)
      shared_peaks <- granges_to_df(gr1[unique(queryHits(overlaps))])
      
      file_path <- sprintf("peak_files/%s_%s_%s_%s_shared_peaks_all.bed", 
                           sample_name, histone_mark, method1, method2)
      write_peaks_to_bed(shared_peaks, file_path)
      cat(sprintf("Written %d shared peaks between %s and %s to %s\n", 
                  nrow(shared_peaks), method1, method2, file_path))
    }
  }
  if (length(methods) >= 3) {
    combinations <- combn(methods, 3, simplify = FALSE)
    
    for (combo in combinations) {
      method1 <- combo[1]
      method2 <- combo[2]
      method3 <- combo[3]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      gr3 <- df_to_granges(method_peaks[[method3]])
      overlaps_12 <- findOverlaps(gr1, gr2)
      shared_gr_12 <- gr1[unique(queryHits(overlaps_12))]
      overlaps_123 <- findOverlaps(shared_gr_12, gr3)
      shared_peaks <- granges_to_df(shared_gr_12[unique(queryHits(overlaps_123))])
      
      file_path <- sprintf("peak_files/%s_%s_%s_%s_%s_shared_peaks_all.bed", 
                           sample_name, histone_mark, method1, method2, method3)
      write_peaks_to_bed(shared_peaks, file_path)
      cat(sprintf("Written %d shared peaks between %s, %s and %s to %s\n", 
                  nrow(shared_peaks), method1, method2, method3, file_path))
    }
  }
  if (length(methods) >= 2) {
    gr_list <- lapply(methods, function(m) df_to_granges(method_peaks[[m]]))
    consensus_gr <- gr_list[[1]]
    for (i in 2:length(gr_list)) {
      overlaps <- findOverlaps(consensus_gr, gr_list[[i]])
      consensus_gr <- consensus_gr[unique(queryHits(overlaps))]
      
      if (length(consensus_gr) == 0) break
    }
    
    if (length(consensus_gr) > 0) {
      consensus_peaks <- granges_to_df(consensus_gr)
      file_path <- sprintf("peak_files/%s_%s_consensus_peaks_all.bed", 
                           sample_name, histone_mark)
      write_peaks_to_bed(consensus_peaks, file_path)
      cat(sprintf("Written %d consensus peaks to %s\n", 
                  nrow(consensus_peaks), file_path))
    }
  }
}
write_results_to_csv <- function(results, method_peak_counts, histone_mark, sample_name, file_path) {
  csv_data <- data.frame(
    Measure = character(),
    Value = integer(),
    stringsAsFactors = FALSE
  )
  for (method in names(method_peak_counts)) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = paste(method, "Total Peaks"),
      Value = method_peak_counts[[method]]
    ))
  }
  unique_measures <- grep("_unique$", names(results), value = TRUE)
  for (measure in unique_measures) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = measure,
      Value = results[[measure]]
    ))
  }
  pairwise_measures <- grep("_shared$", names(results), value = TRUE)
  pairwise_measures <- pairwise_measures[!grepl("_[^_]+_[^_]+_[^_]+_shared$", pairwise_measures)]
  for (measure in pairwise_measures) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = measure,
      Value = results[[measure]]
    ))
  }
  threeway_measures <- grep("_shared$", names(results), value = TRUE)
  threeway_measures <- threeway_measures[grepl("_[^_]+_[^_]+_[^_]+_shared$", threeway_measures)]
  for (measure in threeway_measures) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = measure,
      Value = results[[measure]]
    ))
  }
  if ("consensus" %in% names(results)) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = "Consensus Peaks",
      Value = results[["consensus"]]
    ))
  }
  
  write.csv(csv_data, file = file_path, row.names = FALSE)
  cat(sprintf("Results written to %s\n", file_path))
}

create_histone_upset <- function(results_data, histone_mark, sample_name) {
  results_mark <- results_data$results
  method_peak_counts <- results_data$method_peak_counts
  methods <- names(method_peak_counts)
  total_peaks_per_method <- unlist(method_peak_counts)
  peak_df <- data.frame(matrix(0, nrow = max(sum(total_peaks_per_method), 1), ncol = length(methods)))
  colnames(peak_df) <- methods
  current_row <- 1
  actual_rows_filled <- 0

  if ("consensus" %in% names(results_mark) && results_mark[["consensus"]] > 0) {
    n_consensus <- results_mark[["consensus"]]
    if (current_row + n_consensus - 1 <= nrow(peak_df)) {
      peak_df[current_row:(current_row + n_consensus - 1), ] <- 1
      current_row <- current_row + n_consensus
      actual_rows_filled <- actual_rows_filled + n_consensus
    }
  }
  
  three_way_shared <- grep("_shared$", names(results_mark), value = TRUE)
  three_way_shared <- three_way_shared[grepl("^[^_]+_[^_]+_[^_]+_shared$", three_way_shared)]
  for (key in three_way_shared) {
    methods_shared <- strsplit(sub("_shared$", "", key), "_")[[1]]
    n_shared <- results_mark[[key]]
    if (n_shared > 0 && all(methods_shared %in% methods)) {
      if (current_row + n_shared - 1 <= nrow(peak_df)) {
        peak_df[current_row:(current_row + n_shared - 1), methods_shared] <- 1
        current_row <- current_row + n_shared
        actual_rows_filled <- actual_rows_filled + n_shared
      }
    }
  }
  
  two_way_shared <- grep("_shared$", names(results_mark), value = TRUE)
  two_way_shared <- two_way_shared[!grepl("_[^_]+_[^_]+_[^_]+_shared$", two_way_shared)]
  for (key in two_way_shared) {
    methods_shared <- strsplit(sub("_shared$", "", key), "_")[[1]]
    n_shared <- results_mark[[key]]
    if (n_shared > 0 && all(methods_shared %in% methods)) {
      if (current_row + n_shared - 1 <= nrow(peak_df)) {
        peak_df[current_row:(current_row + n_shared - 1), methods_shared] <- 1
        current_row <- current_row + n_shared
        actual_rows_filled <- actual_rows_filled + n_shared
      }
    }
  }
  
  for (method in methods) {
    unique_key <- paste0(method, "_unique")
    if (unique_key %in% names(results_mark)) {
      n_unique <- results_mark[[unique_key]]
      if (n_unique > 0) {
        if (current_row + n_unique - 1 <= nrow(peak_df)) {
          peak_df[current_row:(current_row + n_unique - 1), method] <- 1
          current_row <- current_row + n_unique
          actual_rows_filled <- actual_rows_filled + n_unique
        }
      }
    }
  }
  
  if (actual_rows_filled > 0) {
    peak_df <- peak_df[1:actual_rows_filled, , drop = FALSE]
  } else {
    cat(sprintf("\nWarning: No overlapping peaks found for %s - %s. Skipping UpSet plot.\n", 
                histone_mark, sample_name))
    return(NULL)
  }
  
  # Create the plot
  tryCatch({
    pdf(sprintf("%s_%s_upset_plot_all.pdf", sample_name, histone_mark), width = 17, height = 10)
    
    print(
      upset(peak_df, methods,
            width_ratio = 0.1,
            min_size = 0,
            n_intersections = 40,
            matrix = (
              intersection_matrix(
                geom = geom_point(size = 3)
              )
              + theme(text = element_text(size = 12))
            ),
            base_annotations = list(
              'Intersection size' = intersection_size(
                text = list(size = 3)
              )
            ),
            set_sizes = (
              upset_set_size(
                geom = geom_bar(fill = "maroon"),
                position = "right"
              )
              + theme(axis.text.x = element_text(angle = 90))
            ),
            themes = upset_modify_themes(
              list(
                intersections_matrix = theme(text = element_text(size = 8)),
                overall = theme(
                  plot.title = element_text(face = "bold", size = 10, hjust = 0.3),
                  text = element_text(size = 10)
                )
              )
            )
      ) +
        ggtitle(sprintf("%s - %s", sample_name, histone_mark)) +
        theme(plot.title = element_text(hjust = 0.3))
    )
    
    dev.off()
    return(peak_df)
  }, error = function(e) {
    cat(sprintf("\nError creating UpSet plot for %s - %s: %s\n", 
                histone_mark, sample_name, e$message))
    return(NULL)
  })
}

run_complete_analysis <- function(file_paths, 
                                  histone_marks = c("H3K27me3", "H3K27ac", "H3K4me3")) {
  all_results <- list()
  upset_results <- list()
  
  all_files <- unlist(file_paths)
  samples <- unique(sapply(all_files, get_sample_name))
  samples <- samples[!is.na(samples)]
  
  for (mark in histone_marks) {
    for (sample in samples) {
      method_count <- 0
      for (method in names(file_paths)) {
        files <- file_paths[[method]]
        mark_sample_files <- files[grep(paste0("H3K", "[^_]+_", sample, "_R"), files)]
        if (length(mark_sample_files) > 0) {
          method_count <- method_count + 1
          cat(sprintf("Found %d files for %s - %s using %s\n", 
                      length(mark_sample_files), sample, mark, method))
        }
      }
      if (method_count < 2) {
        cat(sprintf("Warning: Only %d method(s) found for %s - %s\n", 
                    method_count, sample, mark))
      }
    }
  }
  
  for (mark in histone_marks) {
    all_results[[mark]] <- list()
    upset_results[[mark]] <- list()
    cat(sprintf("\nAnalyzing %s...\n", mark))
    
    for (sample in samples) {
      cat(sprintf("\nProcessing sample: %s\n", sample))
      
      results <- tryCatch({
        analyze_histone_mark_sample(file_paths, mark, sample)
      }, error = function(e) {
        cat(sprintf("Error processing %s - %s: %s\n", mark, sample, e$message))
        NULL
      })
      
      if (!is.null(results)) {
        all_results[[mark]][[sample]] <- results
        upset_results[[mark]][[sample]] <- create_histone_upset(results, mark, sample)
        cat(sprintf("\nResults for %s - %s:\n", sample, mark))
        
        if (!is.null(results$results)) {
          for (method in names(file_paths)) {
            if (paste0(method, "_unique") %in% names(results$results)) {
              cat(sprintf("%s: %d\n", method, results$results[[paste0(method, "_unique")]]))
            }
          }
        
          shared_pairs <- grep("_shared$", names(results$results), value = TRUE)
          shared_pairs <- shared_pairs[!grepl("_[^_]+_[^_]+_[^_]+_shared$", shared_pairs)]
          if (length(shared_pairs) > 0) {
            cat("\nPairwise shared peaks:\n")
            for (key in shared_pairs) {
              methods <- strsplit(sub("_shared$", "", key), "_")[[1]]
              cat(sprintf("%s and %s: %d\n", methods[1], methods[2], results$results[[key]]))
            }
          }
          
          three_way_shared <- grep("_shared$", names(results$results), value = TRUE)
          three_way_shared <- three_way_shared[grepl("^[^_]+_[^_]+_[^_]+_shared$", three_way_shared)]
          if (length(three_way_shared) > 0) {
            cat("\nThree-way shared peaks:\n")
            for (key in three_way_shared) {
              methods <- strsplit(sub("_shared$", "", key), "_")[[1]]
              cat(sprintf("%s, %s and %s: %d\n", methods[1], methods[2], methods[3], results$results[[key]]))
            }
          }
          
          # Print consensus peaks
          if ("consensus" %in% names(results$results)) {
            cat(sprintf("\nConsensus peaks (shared among all methods): %d\n", 
                        results$results[["consensus"]]))
          }
        }
      }
    }
  }
  
  return(list(
    results = all_results,
    upset = upset_results
  ))
}


results <- tryCatch({
  run_complete_analysis(file_paths)
}, error = function(e) {
  cat("Error in analysis:", e$message, "\n")
  NULL
})

