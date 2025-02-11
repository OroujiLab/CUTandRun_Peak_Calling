#!/bin/bash
#SBATCH -J cnr_snr_flex
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -p veryhimem
#SBATCH --time=1-10:00:00
#SBATCH --mem=180G
#SBATCH -e cnr_snr_flex.err
#SBATCH -o cnr_snr_flex.out

# Load required modules
module load bedtools
module load samtools

# Define paths (adjust these according to your file locations)
BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted"
PEAK_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2"
OUTPUT_DIR="${PEAK_DIR}/snr_analysis_all"

mkdir -p $OUTPUT_DIR

# Print initial analysis summary
echo "=== Analysis Summary ===" | tee ${OUTPUT_DIR}/analysis_summary.txt
echo "Date: $(date)" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
echo "" | tee -a ${OUTPUT_DIR}/analysis_summary.txt

# Get unique samples (excluding replicates)
echo "=== Sample Overview ===" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
total_unique_samples=$(find $BAM_DIR -name "*.sorted.bam" | grep -v "IgG" | sed 's/_R[0-9]*\.sorted\.bam$//' | sort -u | wc -l)
echo "Total Unique Samples: $total_unique_samples" | tee -a ${OUTPUT_DIR}/analysis_summary.txt

# List unique samples
echo -e "\nSamples:" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
find $BAM_DIR -name "*.sorted.bam" | grep -v "IgG" | sed 's/_R[0-9]*\.sorted\.bam$//' | sort -u | while read -r line; do
    basename=$(basename "$line")
    echo "- $basename" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
done

# Count marks and their replicates
echo -e "\n=== Marks and Replicates ===" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
find $BAM_DIR -name "*.sorted.bam" | grep -v "IgG" | while read file; do
    basename=$(basename $file .sorted.bam)
    mark=$(echo $basename | cut -d'_' -f2)
    echo "$mark"
done | sort | uniq -c | while read count mark; do
    echo "- $mark: $count replicates" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
done

# Count samples by cell type
echo -e "\n=== Samples by Cell Type ===" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
find $BAM_DIR -name "*.sorted.bam" | grep -v "IgG" | while read file; do
    basename=$(basename $file .sorted.bam)
    celltype=$(echo $basename | cut -d'_' -f3- | sed 's/_R[0-9]*$//')
    echo "$celltype"
done | sort | uniq -c | while read count celltype; do
    echo "- $celltype: $count samples" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
done

echo -e "\n=== Available IgG Controls ===" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
find $BAM_DIR -name "*IgG*.sorted.bam" | while read file; do
    basename=$(basename $file .sorted.bam)
    echo "- $basename" | tee -a ${OUTPUT_DIR}/analysis_summary.txt
done

echo -e "\n=== Starting Analysis ===\n" | tee -a ${OUTPUT_DIR}/analysis_summary.txt

# Function to find matching IgG control
find_igg_control() {
    local basename=$1
    local celltype=$(echo $basename | cut -d'_' -f3-)
    celltype=${celltype%_R*}  # Remove replicate information
    local replicate=$(echo $basename | grep -oP '_R\d+' | sed 's/^_//')
    
    # Look for an exact match first
    local exact_match=$(find $BAM_DIR -name "*IgG*${celltype}*${replicate}*.bam" | head -n 1)
    
    if [ -n "$exact_match" ]; then
        echo $exact_match
    else
        # If no exact match, look for a control with the same cell type but potentially different replicate
        find $BAM_DIR -name "*IgG*${celltype}*.bam" | head -n 1
    fi
}

# Function to get peak file extension
get_peak_extension() {
    local caller=$1
    local mark=$2
    case $caller in
        "GOPEAKS") echo "_peaks.bed" ;;
        "LANCEOTRON") echo "_L-tron.bed" ;;
        "MACS2")
            if [[ $mark == "H3K27me3" ]]; then
                echo "_peaks.broadPeak"
            else
                echo "_peaks.narrowPeak"
            fi
            ;;
        "SEACR") echo "_stringent.stringent.bed" ;;
    esac
}

# Function to standardize peak files to 3 columns
standardize_peak_file() {
    local input_file=$1
    local output_file=$2
    
    # Extract first 3 columns regardless of file format
    cut -f1-3 $input_file > $output_file
    
    # Check if file is empty or has incorrect format
    if [[ ! -s $output_file ]]; then
        echo "Warning: Standardized peak file is empty: $output_file"
        return 1
    fi
    
    # Count peaks in standardized file
    local peak_count=$(wc -l < $output_file)
    echo "Standardized peaks: $peak_count"
}

# Function to calculate SNR with detailed statistics
calculate_snr() {
    local cnr_bam=$1
    local igg_bam=$2
    local peaks=$3
    local output_prefix=$4
    local mark=$5
    local caller=$6
    local sample_name=$7
    
    # Standardize peak file
    local std_peaks="${output_prefix}.standardized_peaks_all.bed"
    standardize_peak_file $peaks $std_peaks
    
    # Create background regions
    samtools view -H $cnr_bam | grep "@SQ" | cut -f 2,3 | sed 's/SN://g' | sed 's/LN://g' > ${output_prefix}.genome
    
    # Generate random regions of same number and size as peaks
    bedtools shuffle -i $std_peaks -g ${output_prefix}.genome > ${output_prefix}.background_all.bed
    
    # Calculate coverage for cnr and IgG in peak regions
    bedtools coverage -a $std_peaks -b $cnr_bam | sort -k4,4n > ${output_prefix}.cnr_peak_coverage.sorted_all.txt
    bedtools coverage -a $std_peaks -b $igg_bam | sort -k4,4n > ${output_prefix}.igg_peak_coverage.sorted_all.txt
    
    # Calculate coverage for cnr and IgG in background regions
    bedtools coverage -a ${output_prefix}.background_all.bed -b $cnr_bam | sort -k4,4n > ${output_prefix}.cnr_background_coverage.sorted_all.txt
    bedtools coverage -a ${output_prefix}.background_all.bed -b $igg_bam | sort -k4,4n > ${output_prefix}.igg_background_coverage.sorted_all.txt
    
    # Calculate statistics for each dataset
    for dataset in cnr_peak igg_peak cnr_background igg_background; do
        awk -v mark="$mark" -v caller="$caller" -v dataset="$dataset" -v sample="$sample_name" '
        BEGIN {
            sum=0; count=0; min=1e30; max=-1e30;
        }
        {
            sum+=$4;
            count++;
            if($4 < min) min=$4;
            if($4 > max) max=$4;
            values[count]=$4;
        }
        END {
            if(count > 0) {
                mean=sum/count;
                
                # Calculate median
                if(count%2) {
                    median=values[int(count/2)+1];
                } else {
                    median=(values[count/2]+values[count/2+1])/2;
                }
                
                print mark "\t" caller "\t" sample "\t" dataset "_Count\t" count;
                print mark "\t" caller "\t" sample "\t" dataset "_Mean\t" mean;
                print mark "\t" caller "\t" sample "\t" dataset "_Median\t" median;
                print mark "\t" caller "\t" sample "\t" dataset "_Min\t" min;
                print mark "\t" caller "\t" sample "\t" dataset "_Max\t" max;
            }
        }' ${output_prefix}.${dataset}_coverage.sorted_all.txt > ${output_prefix}.${dataset}_stats_all.txt
    done
    
    # Combine all stats files
    cat ${output_prefix}.*_stats_all.txt > ${output_prefix}.stats_all.txt
    
    # Calculate final SNR values
    awk -v mark="$mark" -v caller="$caller" -v sample="$sample_name" '
    /cnr_peak_Mean/ {cnr_peak_mean=$5}
    /igg_peak_Mean/ {igg_peak_mean=$5}
    /cnr_background_Mean/ {cnr_bg_mean=$5}
    /igg_background_Mean/ {igg_bg_mean=$5}
    /cnr_peak_Median/ {cnr_peak_median=$5}
    /igg_peak_Median/ {igg_peak_median=$5}
    /cnr_background_Median/ {cnr_bg_median=$5}
    /igg_background_Median/ {
        igg_bg_median=$5;
        
        # Calculate SNR values with IgG correction
        if(igg_peak_mean > 0 && igg_bg_mean > 0) {
            cnr_igg_mean_snr=(cnr_peak_mean/igg_peak_mean)/(cnr_bg_mean/igg_bg_mean);
        } else {
            cnr_igg_mean_snr=0;
        }
        
        if(igg_peak_median > 0 && igg_bg_median > 0) {
            cnr_igg_median_snr=(cnr_peak_median/igg_peak_median)/(cnr_bg_median/igg_bg_median);
        } else {
            cnr_igg_median_snr=0;
        }
        
        print mark "\t" caller "\t" sample "\tcnr_IgG_Mean_SNR\t" cnr_igg_mean_snr;
        print mark "\t" caller "\t" sample "\tcnr_IgG_Median_SNR\t" cnr_igg_median_snr;
    }' ${output_prefix}.stats_all.txt >> ${output_prefix}.stats_all.txt
    
    # Save top 10 highest coverage peaks
    head -n 10 ${output_prefix}.cnr_peak_coverage.sorted_all.txt > ${output_prefix}.top10_peaks_all.txt
}

# Get total number of samples for progress tracking
total_samples=$(find $BAM_DIR -name "*.sorted.bam" | grep -v "IgG" | wc -l)
current_sample=0

# Find all cnr BAM files (excluding IgG controls)
find $BAM_DIR -name "*.sorted.bam" | grep -v "IgG" | while read cnr_bam; do
    basename=$(basename $cnr_bam .sorted.bam)
    ((current_sample++))
    
    # Print processing status
    echo -e "\nProcessing Sample $current_sample of $total_samples: $basename"
    
    # Extract mark and cell type information
    mark=$(echo $basename | cut -d'_' -f2)
    celltype=$(echo $basename | cut -d'_' -f3- | sed 's/_R[0-9]*$//')
    replicate=$(echo $basename | grep -oP '_R\d+' | sed 's/^_//')
    
    echo "Mark: $mark"
    echo "Cell Type: $celltype"
    echo "Replicate: $replicate"
    
    # Find matching IgG control
    igg_bam=$(find_igg_control $basename)
    
    if [ -z "$igg_bam" ]; then
        echo "Warning: No matching IgG control found for $basename - using CNR BAM as IgG"
        igg_bam=$cnr_bam
    else
        echo "Using IgG control: $(basename $igg_bam)"
    fi
    
    # Process with each peak caller
    for CALLER in "GOPEAKS" "LANCEOTRON" "MACS2" "SEACR"; do
        PEAK_EXT=$(get_peak_extension $CALLER $mark)
        PEAK_FILE="${PEAK_DIR}/${CALLER}/${basename}${PEAK_EXT}"
        
        if [ ! -f "$PEAK_FILE" ]; then
            echo "Warning: Peak file not found: $PEAK_FILE"
            continue
        fi
        
        OUTPUT_PREFIX="${OUTPUT_DIR}/${basename}_${CALLER}"
        
        echo "Processing with ${CALLER}"
        calculate_snr "$cnr_bam" "$igg_bam" "$PEAK_FILE" "$OUTPUT_PREFIX" "$mark" "$CALLER" "$basename"
    done
done

# Combine all results
echo -e "Mark\tCaller\tSample\tMetric\tValue" > ${OUTPUT_DIR}/snr_summary_all.tsv
cat ${OUTPUT_DIR}/*stats_all.txt >> ${OUTPUT_DIR}/snr_summary_all.tsv

# Create focused SNR summary
echo -e "Mark\tCaller\tSample\tSNR_Type\tValue" > ${OUTPUT_DIR}/snr_comparison_all.tsv
awk -F'\t' '
    /_SNR/ {
        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5
    }
' ${OUTPUT_DIR}/snr_summary_all.tsv | sort -k1,1 -k2,2 -k3,3 >> ${OUTPUT_DIR}/snr_comparison_all.tsv

echo -e "\nAnalysis completed"
echo "Summary files:"
echo "1. ${OUTPUT_DIR}/analysis_summary.txt"
echo "2. ${OUTPUT_DIR}/snr_summary_all.tsv"
echo "3. ${OUTPUT_DIR}/snr_comparison_all.tsv"
