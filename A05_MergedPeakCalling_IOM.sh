#!/bin/bash

# Purpose:
#    Run peak calling on merged fragments
#    Perform IOM on Peaks


# 1. Peak calling with MACS3
### PC_IOM.sh
# (the .bed.gz is the fragment file exported in the previous step)

set -e

EXP_name="${1}";
IOM_threads="${2}";

echo "Peak Calling on ${EXP_name}"

macs3 callpeak --treatment ${EXP_name}.bed.gz \
--format BEDPE \
--gsize mm \
--nomodel \
--nolambda \
--keep-dup all \
--qvalue 0.01 \
--call-summits \
--outdir ${EXP_name}_macs3_Q01 \
--name ${EXP_name} \
--verbose 2

#macs3 hmmratac --input ${EXP_name}.bed.gz \
#--format BEDPE \
#--outdir macs3_hmmratac \
#--name ${EXP_name} \

echo "Iterative Overlap Merging on ${EXP_name}_macs3_Q01/${EXP_name}_summits.bed"

Rscript IOM.R "$(pwd)" 500 "${EXP_name}_macs3_Q01/${EXP_name}_summits.bed" "${EXP_name}_ITMPeaks.bed" ${IOM_threads}
gzip -c ${EXP_name}_ITMPeaks.bed > ${EXP_name}_ITMPeaks.bed.gz

echo "Peak Calling and IOM done for ${EXP_name}"

#ENV:
#mm env create --file PCIOM.yml
#or
#mm create -n PCIOM python=3.11
#mmac PCIOM
#pip install MACS3
#mm install r-base r-essentials r-data.table bioconductor-biocinstaller bioconductor-genomicranges

#Activate mPCIOM env before running
# bash A05_MergedPeakCalling_IOM.sh