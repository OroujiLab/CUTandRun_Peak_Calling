![](images/Readme_2.png)

# Benchmarking Peak Calling Methods for CUT&RUN:

Cleavage Under Targets and Release Using Nuclease (CUT&RUN) has rapidly gained prominence as an effective approach for mapping protein-DNA interactions, especially histone modifications, offering substantial improvements over conventional chromatin immunoprecipitation sequencing (ChIP-seq). However, the effectiveness of this technique is contingent upon accurate peak identification, necessitating the use of optimal peak calling methods tailored to the unique characteristics of CUT&RUN data. Here, we benchmark four prominent peak calling tools - MACS2, SEACR, GoPeaks, and LanceOtron - evaluating their performance in identifying peaks from CUT&RUN datasets. Our analysis utilizes in-house data of three histone marks (H3K4me3, H3K27ac, and H3K27me3) from mouse brain tissue, as well as samples from the 4DNucleome database. We systematically assess these tools based on parameters such as the number of peaks called, peak length distribution, signal enrichment, computational efficiency, and reproducibility across biological replicates. Our findings reveal substantial variability in peak calling efficacy, with each method demonstrating distinct strengths in sensitivity, precision, and applicability depending on the histone mark in question. These insights provide a comprehensive evaluation that will assist in selecting the most suitable peak caller for high-confidence identification of regions of interest in CUT&RUN experiments, ultimately enhancing the study of chromatin dynamics and transcriptional regulation.


# Citation:

If you use the scripts or methods from this repository, please cite the preprint:

Benchmarking Peak Calling Methods for CUT&RUN

Amin Nooranikhojasteh, Ghazaleh Tavallaee, and Elias Orouji 
bioRxiv, 2024. DOI: https://doi.org/10.1101/2024.11.13.622880


