#!/bin/bash

# Purpose:
#    Run peak calling after 10x cellcalling and QC
#    Perform IOM on Peaks

# 1. Peak calling with MACS3

# Use a for loop to iterate over files ending in "bed.gz"
# (the .bed.gz are fragments exported in the previous step)
for file in ./*.bed.gz; do
    # Get the filename without extension
    filename=$(basename "$file" .bed.gz)
    echo "Peak calling: $filename"
    
    macs3 callpeak --treatment $file \
        --format BEDPE \
        --gsize mm \
        --nomodel \
        --nolambda \
        --keep-dup all \
        --qvalue 0.01 \
        --call-summits \
        --outdir ${filename}_macs3_Q01 \
        --name $filename \
        --verbose 2
    echo "Peak calling done: $filename"
done
###

# 2. Iterative Overlap Merging of the Peaks

# Iterate through folders ending with "_macs3_Q01"
# and perform IterativeOverlapMerging of the peaks
# https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html

for folder in *_macs3_Q01; do
    if [ -d "$folder" ]; then
        # Get folder name without "_macs3_Q01"
        folder_name="${folder%_macs3_Q01}"
        echo "Peak Merging: $folder_name"

        # Change directory to the current folder
        cp IOM.R "$folder"
        cd "$folder"

        # Launch the R script with the specified arguments
        Rscript IOM.R "$(pwd)" 500 "${folder_name}_summits.bed" "${folder_name}_ITMPeaks.bed 20"
        echo "Peak Merging done: $folder_name"

        # Move back to the parent directory
        rm IOM.R
        cd ..
    fi
done
###

#ENV:
#mm env create --file PCIOM.yml
#or
#mm create -n PCIOM python=3.11
#mmac PCIOM
#pip install MACS3
#mm install r-base r-essentials r-data.table bioconductor-biocinstaller bioconductor-genomicranges

#Activate mPCIOM env before running
# bash A02_BatchPeakCalling_IOM.sh
