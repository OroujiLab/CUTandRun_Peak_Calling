#!/bin/bash


BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted"
OUTPUT_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/LANCEOTRON"
BIGWIG_DIR="${OUTPUT_DIR}/bigwig"

mkdir -p "$OUTPUT_DIR"

run_lanceotron() {
    local bigwig_file=$1
    local output_prefix=$2
    local sample_dir="${OUTPUT_DIR}/${output_prefix}"

    if [ -d "$sample_dir" ] && [ -f "${sample_dir}/${output_prefix}_L-tron.bed" ]; then
        echo "Output for $bigwig_file already exists in ${sample_dir}. Skipping..."
        return
    fi

    mkdir -p "$sample_dir"
    
    echo "Processing: $bigwig_file with LANCEOTRON"    
    lanceotron callPeaks "$bigwig_file" -f "${sample_dir}/${output_prefix}_L-tron"
    
    echo "LANCEOTRON completed for ${output_prefix}"
}

for bam_file in ${BAM_DIR}/*.bam; do
    if [[ $bam_file == *.bai ]] || [[ $bam_file == *IgG* ]]; then
        continue
    fi
    
    base_name=$(basename "$bam_file" .sorted.bam)   
    if [ -f "${BIGWIG_DIR}/${base_name}.bw" ]; then
        echo "Found existing bigWig file for: $bam_file"
        
        run_lanceotron "${BIGWIG_DIR}/${base_name}.bw" "$base_name"
    else
        echo "No corresponding bigWig file found for: $bam_file. Skipping..."
    fi
done

echo "All available bigWig files have been processed with LANCEOTRON."
