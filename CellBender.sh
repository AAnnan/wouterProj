#!/bin/bash

Dir10x='/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
PathMat='/outs/raw_feature_bc_matrix.h5'
n_threads='12'

# Get a list of samples, they start with 'WK'
Samples=($(ls "$Dir10x" | grep '^WK'))
#Samples=("WK-1585_Regen_Day3_AP_BL6" "WK-1585_INTACT_AP_BL6_Contrl")

for sample in "${Samples[@]}"; do
    n_cell=$(awk -F',' '{sum += $4} END {print sum}' "$Dir10x$sample/outs/per_barcode_metrics.csv")
    input_mat="$Dir10x$sample/$PathMat"
    output_h5="${sample}_CellBender.h5"

    cellbender remove-background \
        --input "$input_mat" \
        --output "$output_h5" \
        --expected-cells "$n_cell" \
        --cpu-threads "$n_threads" \
        --checkpoint-mins 180

done
