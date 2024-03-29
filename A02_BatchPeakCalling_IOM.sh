#!/bin/bash

set -e

# Purpose:
#    Run peak calling after 10x cellcalling and QC
#    Perform IOM on Peaks

# Usage:
#ENV:
#mm env create --file PCIOM.yml
#mmac PCIOM
#pip install MACS3 numpy pandas scipy statsmodels
#or
#mm create -n PCIOM python=3.11
#mmac PCIOM
#pip install MACS3 numpy pandas scipy statsmodels
#mm install r-base r-essentials r-data.table bioconductor-biocinstaller bioconductor-genomicranges

#Activate PCIOM env before running bash A02_BatchPeakCalling_IOM.sh
# Have IOM.R in folder

#Amulet (https://github.com/UcarLab/AMULET)
# Changes: 
#                   FragmentFileOverlapCounter.py line 231  'is__cell_barcode' -> 'is_cell'
#                   FragmentFileOverlapCounter.py line 193  np.object -> object
#                   peakoverlap.py line 193 np.object -> object
#                   AMULET.py line 102,133,138  np.object -> object


# Script

# 1. Peak calling with MACS3

# Use a for loop to iterate over files ending in "tsv.gz"
# (the .tsv.gz are fragments exported in the previous step)
for file in ./*.tsv.gz; do
    # Get the filename without extension
    filename=$(basename "$file" .tsv.gz)
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
        echo Launching IOM.R "$(pwd)" 500 "${folder_name}_summits.bed" "${folder_name}_ITMPeaks.bed" 20
        Rscript IOM.R "$(pwd)" 500 "${folder_name}_summits.bed" "${folder_name}_ITMPeaks.bed" 20
        echo "Peak Merging done: $folder_name"

        # Move back to the parent directory
        rm IOM.R
        cd ..
    fi
done
###


# 3. AMULET count based method for detecting multiplets from snATAC-seq data.

# Use a for loop to iterate over files ending in "tsv.gz"
# (the .tsv.gz are fragments exported in the previous step)
for file in ./*.tsv.gz; do
    # Get the filename without extension
    filename=$(basename "$file" .tsv.gz)
    echo "Detecting multiplets: $filename"
    
    mkdir -p ${filename}_amulet
    
    ./amulet/AMULET.sh $file /mnt/ndata/daniele/wouter/Processed/CellRangerArc/${filename}/outs/per_barcode_metrics.csv ./amulet/mouse_chromosomes_noxy.txt ./amulet/mm10-blacklist.v2.bed ${filename}_amulet ./amulet/.

    echo -e "DONE: $filename\n"
done
###

