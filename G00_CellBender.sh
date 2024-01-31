#!/bin/bash

##############################
## SCRIPT RUNNING ON UNIL GPUs
##############################
#! /bin/bash

## Resource Allocation
#SBATCH --time=3-00:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

## job metadata
#SBATCH --job-name="CellBender"
#SBATCH --mail-user=ahrmad.annan@unil.ch
#SBATCH --mail-type=end,fail
set -e

source ~/micromamba/etc/profile.d/micromamba.sh

micromamba activate cellbender

Dir10x='/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
PathRawMat='/outs/raw_feature_bc_matrix.h5'
Samples=("WK-1350_I-1_AP" "WK-1350_I-1_LP" "WK-1350_I-2_AP" "WK-1350_I-2_DP" "WK-1350_R3-1_AP" "WK-1350_R3-1_DP" "WK-1350_R3-1_LP" "WK-1350_R3-1_VP" "WK-1350_R3-2_AP" "WK-1350_R3-2_DP" "WK-1350_R3-2_LP" "WK-1350_R3-2_VP" "WK-1384_BL6_Intact_AP_2_SLT" "WK-1501_BL6_INTACT_AP_Test3" "WK-1501_BL6_INTACT_AP_Test3_SORT" "WK-1580_BL6_AP_RegenDay1" "WK-1580_BL6_AP_RegenDay2" "WK-1585_Castrate_Day28_AP_BL6" "WK-1585_INTACT_AP_BL6_Citrate" "WK-1585_INTACT_AP_BL6_Contrl" "WK-1585_Regen_Day3_AP_BL6")
#mapfile -t Samples < "WK_samples"
#Samples=("WK-1501_BL6_INTACT_AP_Test3_SORT" "WK-1580_BL6_AP_RegenDay1")

mkdir raws
for sample in "${Samples[@]}"; do
    scp ahrmad@electron.unil.ch:/mnt/ndata/daniele/wouter/Processed/CellRangerArc/${sample}/outs/raw_feature_bc_matrix.h5 raws/${sample}_raw_feature_bc_matrix.h5
    scp ahrmad@electron.unil.ch:/mnt/ndata/daniele/wouter/Processed/CellRangerArc/${sample}/outs/per_barcode_metrics.csv raws/${sample}_per_barcode_metrics.csv

    n_cell=$(awk -F',' '{sum += $4} END {print sum}' "raws/${sample}_per_barcode_metrics.csv")
    input_mat="raws/${sample}_raw_feature_bc_matrix.h5"
    output_h5="${sample}_CellBender.h5"

    cellbender remove-background \
        --input "$input_mat" \
        --output "$output_h5" \
        --expected-cells "$n_cell" \
        --checkpoint-mins 1440 \
        --cuda

done
micromamba deactivate 

#sbatch CB.sh

##############################
## FILE CHECK
##############################

#!/bin/bash

# Your array
Dir10x='/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
PathRawMat='/outs/raw_feature_bc_matrix.h5'
mapfile -t Samples < "WK_samples"
elements_with_missing_files=()

# Function to check file existence
check_files_exist() {
    local element="$1"
    local extensions=("_CellBender.h5" "_CellBender.log" "_CellBender.pdf" "_CellBender_cell_barcodes.csv" "_CellBender_filtered.h5" "_CellBender_metrics.csv" "_CellBender_posterior.h5" "_CellBender_report.html")
    local missing_files=()

    for ext in "${extensions[@]}"; do
        file="${element}${ext}"
        if [ ! -e "$file" ]; then
            missing_files+=("$file")
        fi
    done

    if [ ${#missing_files[@]} -eq 0 ]; then
        echo "All files exist for element '$element'"
    else
        echo "Missing files for element '$element': ${missing_files[@]}"
        elements_with_missing_files+=("$element")
    fi
}

# Iterate through the array and check files
for element in "${Samples[@]}"; do
    check_files_exist "$element"
done

# Output elements with missing files
if [ ${#elements_with_missing_files[@]} -gt 0 ]; then
    echo "Elements with one or more missing files:"
    echo "${elements_with_missing_files[@]}"
else
    echo "All elements have all required files."
fi

##############################
### ORIGINAL CPU SCRIPT
##############################

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





