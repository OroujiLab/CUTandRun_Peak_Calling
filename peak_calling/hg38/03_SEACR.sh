#!/bin/bash


BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted"
OUTPUT_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/SEACR"
BEDGRAPH_DIR="${OUTPUT_DIR}/bedgraphs"

mkdir -p "$OUTPUT_DIR" "$BEDGRAPH_DIR"

find_igg_control() {
    local basename=$1
    local celltype=$(echo $basename | cut -d'_' -f3-)
    celltype=${celltype%_R*}  # Remove replicate information
    local replicate=$(echo $basename | grep -oP '_R\d+' | sed 's/^_//')
    
    
    local exact_match=$(find $BAM_DIR -name "*IgG*${celltype}*${replicate}*.bam" | head -n 1)
    
    if [ -n "$exact_match" ]; then
        echo $exact_match
    else
        # If no exact match, look for a control with the same cell type but potentially different replicate
        find $BAM_DIR -name "*IgG*${celltype}*.bam" | head -n 1
    fi
}

run_seacr() {
    local bam_file=$1
    local output_prefix=$2
    local control_file=$3

    echo "Processing: $bam_file with SEACR"
    
    if [ -f "${OUTPUT_DIR}/${output_prefix}_relaxed.peaks" ] && [ -f "${OUTPUT_DIR}/${output_prefix}_stringent.peaks" ]; then
        echo "Output for $bam_file already exists. Skipping..."
        return
    fi
    
    # Convert BAM to bedgraph because SEACR uses that for peak calling
    bedtools genomecov -bg -ibam "$bam_file" > "${BEDGRAPH_DIR}/${output_prefix}.bedgraph"
    
    if [ -n "$control_file" ]; then
        bedtools genomecov -bg -ibam "$control_file" > "${BEDGRAPH_DIR}/${output_prefix}_control.bedgraph"
        bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_prefix}.bedgraph" \
          "${BEDGRAPH_DIR}/${output_prefix}_control.bedgraph" \
          norm relaxed "${OUTPUT_DIR}/${output_prefix}_relaxed"
        bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_prefix}.bedgraph" \
          "${BEDGRAPH_DIR}/${output_prefix}_control.bedgraph" \
          norm stringent "${OUTPUT_DIR}/${output_prefix}_stringent"
    else
        bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_prefix}.bedgraph" \
          0.01 non relaxed "${OUTPUT_DIR}/${output_prefix}_relaxed"
        bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_prefix}.bedgraph" \
          0.01 non stringent "${OUTPUT_DIR}/${output_prefix}_stringent"
    fi
    
    echo "SEACR completed for ${output_prefix}"
}


for bam_file in ${BAM_DIR}/*.bam; do
    if [[ $bam_file == *.bai ]] || [[ $bam_file == *IgG* ]]; then
        continue
    fi
    
    base_name=$(basename "$bam_file" .sorted.bam)    
    control_file=$(find_igg_control "$base_name")
    
    echo "Processing file: $bam_file"
    if [ -n "$control_file" ]; then
        echo "Using control: $control_file"
    else
        echo "No matching control found"
    fi
    
    run_seacr "$bam_file" "$base_name" "$control_file"
done

echo "All BAM files have been processed with SEACR."
