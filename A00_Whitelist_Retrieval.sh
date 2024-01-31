#!/bin/bash

set -e

##Whitelists retrieval

# Purpose:
#    To get the barcodes of cells called by 10x pipeline
#    Used as a whitelist filter in the import by SNAPATAC2

# Variables
base_directory="${1}";
output_directory="${2}";

if [ ${#@} -lt 2 ] ; then
    printf '\nUsage:\n';
    printf '    A00_Data_Retrieval.sh \\\n';
    printf '        base10x_directory \\\n';
    printf '        output_directory \\\n';
    exit 1
fi
#base_directory="/mnt/ndata/daniele/wouter/Processed/CellRangerArc"
#output_directory="/mnt/etemp/ahrmad/wouter/batch_ATAC"

# Download Blacklist
wget -N https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz

# Loop over directories
for directory in "$base_directory"/WK*; do
   echo "$directory"
   if [ -d "$directory" ]; then
       csv_file=("$directory"/outs/per_barcode_metrics.csv)
       summ_file=("$directory"/outs/summary.csv)
       

       mkdir -p "$output_directory/whitelists"
       wl="$output_directory/whitelists/$(basename "$directory")"
       wlLOG="$output_directory/whitelists/$(basename "$directory")_LOG"
       awk -F ',' '$4 == 1 {print $1}' "$csv_file" > "$wl"

       echo $(basename "$directory")",Total Barcodes,Cell Barcodes,ATAC Sequenced read pairs,GEX Sequenced read pairs" >> "${wlLOG}"
       echo 0,$(wc -l "$csv_file" | cut -f1 -d" "),$(wc -l "$wl" | cut -f1 -d" "),$(cut -f23,45 -d',' "$summ_file"| sed -n '2 p') >> "${wlLOG}"

       echo $(basename "$directory") whitelist done
   fi
done

