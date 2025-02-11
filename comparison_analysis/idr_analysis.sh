#!/bin/bash

BASE_DIR="/Aminnn/CNR/EpigenomeLab/Four_DN/results"
SNR_DIR="${BASE_DIR}/snr_analysis"
OUTPUT_DIR="${BASE_DIR}/idr_analysis"

mkdir -p $OUTPUT_DIR

process_peaks() {
    local mark=$1
    local celltype=$2
    local caller=$3
    echo "Processing $mark $celltype $caller"    
    mapfile -t files < <(find $SNR_DIR -name "*_${mark}_${celltype}_R*_${caller}.chip_peak_coverage.sorted_all.txt" | sort -V)
    
    if [ ${#files[@]} -lt 2 ]; then
        echo "Skipping $mark $celltype $caller - fewer than 2 replicates found"
        return
    fi
    
    echo "Found ${#files[@]} replicates"
    tmp_dir="${OUTPUT_DIR}/tmp_${mark}_${celltype}_${caller}"
    mkdir -p $tmp_dir
    for file in "${files[@]}"; do
        rep_num=$(echo "$file" | grep -o "_R[0-9]*_" | sed 's/_R\([0-9]*\)_/\1/')
        echo "Processing replicate $rep_num from file: $file"
        
        awk -v OFS="\t" '{print $1,$2,$3,"peak_"NR,$4,".",$4,1,1,0}' "$file" > "${tmp_dir}/rep${rep_num}.narrowPeak"
    done
    
    for ((i=0; i<${#files[@]}-1; i++)); do
        j=$((i+1))
        rep1_num=$(echo "${files[$i]}" | grep -o "_R[0-9]*_" | sed 's/_R\([0-9]*\)_/\1/')
        rep2_num=$(echo "${files[$j]}" | grep -o "_R[0-9]*_" | sed 's/_R\([0-9]*\)_/\1/')
        
        output="${OUTPUT_DIR}/${mark}_${celltype}_${caller}_R${rep1_num}_vs_R${rep2_num}.idr"
        
        echo "Running IDR on R${rep1_num} vs R${rep2_num}"
        idr --samples "${tmp_dir}/rep${rep1_num}.narrowPeak" "${tmp_dir}/rep${rep2_num}.narrowPeak" \
            --input-file-type narrowPeak \
            --rank signal.value \
            --output-file "$output" \
            --plot \
            --log-output-file "${output}.log"
        
        
        if [ -f "idr.png" ]; then
            mv idr.png "${OUTPUT_DIR}/${mark}_${celltype}_${caller}_R${rep1_num}_vs_R${rep2_num}.png"
        fi
        if [ -f "idr.pdf" ]; then
            mv idr.pdf "${OUTPUT_DIR}/${mark}_${celltype}_${caller}_R${rep1_num}_vs_R${rep2_num}.pdf"
        fi
    done
}


find $SNR_DIR -name "*chip_peak_coverage.sorted_all.txt" | while read file; do
    basename=$(basename "$file" .chip_peak_coverage.sorted_all.txt)
    mark=$(echo "$basename" | cut -d'_' -f2)
    celltype=$(echo "$basename" | cut -d'_' -f3)
    caller=$(echo "$basename" | rev | cut -d'_' -f1 | rev)
    key="${mark}_${celltype}_${caller}"
    
    
    if [[ ! -f "${OUTPUT_DIR}/.processed_${key}" ]]; then
        process_peaks "$mark" "$celltype" "$caller"
        touch "${OUTPUT_DIR}/.processed_${key}"
    fi
done

rm -f ${OUTPUT_DIR}/.processed_*
#echo "IDR is completed"
