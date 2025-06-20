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
    id <- parts[1]  # Gets "4DN..."
    histone_idx <- grep("H3K", parts)[1]
    replicate_idx <- grep("^R\\d+", parts)[1]
    
    if (is.na(histone_idx) || is.na(replicate_idx)) {
      stop("Invalid filename format")
    }
    
    histone <- parts[histone_idx]  # Gets "histone mark"
    sample_name <- paste(parts[(histone_idx + 1):(replicate_idx - 1)], collapse = "_")  # Gets "H1_ESC"
    replicate <- parts[replicate_idx]  # Gets "replicate number"
    
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
    warning(sprintf("impossible to get sample name from %s: %s", filename, e$message))
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


create_venn_diagram <- function(method_peaks, histone_mark, sample_name, output_dir) {
  dir.create(file.path(output_dir, sample_name, histone_mark), recursive = TRUE, showWarnings = FALSE)
  peak_granges <- lapply(method_peaks, function(peaks) {
    GRanges(seqnames = peaks$chr,
            ranges = IRanges(start = peaks$start, end = peaks$end))
  })
  
  all_regions <- GenomicRanges::reduce(unlist(GRangesList(peak_granges)))
  
  venn_lists <- list()
  overlap_counts <- list()
  
  for(method in names(peak_granges)) {
    if(!is.null(peak_granges[[method]])) {
      overlaps <- findOverlaps(peak_granges[[method]], all_regions)
      venn_lists[[method]] <- unique(subjectHits(overlaps))
      
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
  
  cat("\nOverlap Statistics for", histone_mark, "-", sample_name, ":\n")
  for(comparison in names(overlap_counts)) {
    cat(sprintf("%s: %d peaks\n", comparison, overlap_counts[[comparison]]))
  }
  
  venn_plot <- ggVennDiagram(venn_lists) +
    scale_fill_viridis() +
    labs(title = paste(histone_mark, "-", sample_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  pdf_file <- file.path(output_dir, sample_name, histone_mark, "venn_diagram.pdf")
  pdf(pdf_file, width = 10, height = 8)
  print(venn_plot)
  dev.off()
  
  return(venn_plot)
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
  
  create_venn_diagram(method_peaks, histone_mark, sample_name, output_dir)
}


generate_peak_analyses <- function(file_paths, output_dir = "peak_analyses") {
  histone_marks <- c("H3K27ac", "H3K27me3", "H3K4me3")
  samples <- unique(sapply(unlist(file_paths), get_sample_name))
  samples <- samples[!is.na(samples)]
  
  cat("Processing the following samples:", paste(samples, collapse=", "), "\n")
  cat("For histone marks:", paste(histone_marks, collapse=", "), "\n\n")
  
  dir.create(output_dir, showWarnings = FALSE)
  
  all_plots <- list()
  
  for (sample in samples) {
    cat(sprintf("\nProcessing sample: %s\n", sample))
    for (mark in histone_marks) {
      cat(sprintf("\nProcessing mark: %s\n", mark))
      plot <- analyze_histone_mark_sample(file_paths, mark, sample, output_dir)
      if (!is.null(plot)) {
        all_plots[[paste(sample, mark, sep="_")]] <- plot
      }
    }
  }
  
  if (length(all_plots) > 0) {
    pdf(file.path(output_dir, "all_venn_diagrams.pdf"), width = 10, height = 8)
    for (plot in all_plots) {
      print(plot)
    }
    dev.off()
  }
  
  #cat("\nAnalyses completed. Results saved in:", output_dir, "\n")
}

# Run 
generate_peak_analyses(file_paths)
