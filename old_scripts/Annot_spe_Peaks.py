import matplotlib.pyplot as plt
import snapatac2 as snap
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
from matplotlib_venn import venn3
import os
import csv

Experiment='Wouter21_SING_CB'
lum_ext = '_annot_Lum_newSigs.h5ad'

adata = sc.read_h5ad(Experiment + lum_ext)
adata_L1 = adata[np.isin(adata.obs['Annotation'],['L1'])].copy()
adata_L2 = adata[np.isin(adata.obs['Annotation'],['L2'])].copy()

adata_CD28 = snap.read('../../Wouter21_SING/WK-1585_Castrate_Day28_AP_BL6_qcTOREMOVE.h5ad').to_memory()

adata_L1_CD28 = adata_L1[np.isin(adata_L1.obs['timePoint'],'Day28')].copy()
adata_CD28_L1 = adata_CD28[np.isin(adata_CD28.obs.index,adata_L1_CD28.obs.index)].copy()
adata_CD28_L1.obs['TYPE'] = ['CD28_L1']*adata_CD28_L1.n_obs
snap.ex.export_fragments(adata_CD28_L1, groupby='TYPE', prefix='CD28_L1_FRAGS', suffix='.tsv.gz')

adata_L2_CD28 = adata_L2[np.isin(adata_L2.obs['timePoint'],'Day28')].copy()
adata_CD28_L2 = adata_CD28[np.isin(adata_CD28.obs.index,adata_L2_CD28.obs.index)].copy()
adata_CD28_L2.obs['TYPE'] = ['CD28_L2']*adata_CD28_L2.n_obs
snap.ex.export_fragments(adata_CD28_L2, groupby='TYPE', prefix='CD28_L2_FRAGS', suffix='.tsv.gz')

adata_L1_I = adata_L1[np.isin(adata_L1.obs['timePoint'],'Intact')].copy()
adata_L2_I = adata_L2[np.isin(adata_L2.obs['timePoint'],'Intact')].copy()

for a,b in [('WK-1501_BL6_INTACT_AP_Test3','WK-1501-Test3'), ('WK-1585_INTACT_AP_BL6_Citrate','WK-1585-Citrate'), ('WK-1585_INTACT_AP_BL6_Contrl','WK-1585-Contrl')]:
	adata_I = snap.read(f'../../Wouter21_SING/{a}_qcTOREMOVE.h5ad').to_memory()
	adata_I_L1 = adata_I[np.isin(adata_I.obs.index,adata_L1_I.obs.index)].copy()
	adata_I_L1.obs['TYPE'] = ['I_L1']*adata_I_L1.n_obs
	snap.ex.export_fragments(adata_I_L1, groupby='TYPE', prefix=f'{b}_I_L1_FRAGS', suffix='.tsv.gz')

	adata_I_L2 = adata_I[np.isin(adata_I.obs.index,adata_L2_I.obs.index)].copy()
	adata_I_L2.obs['TYPE'] = ['I_L2']*adata_I_L2.n_obs
	snap.ex.export_fragments(adata_I_L2, groupby='TYPE', prefix=f'{b}_I_L2_FRAGS', suffix='.tsv.gz')



zcat *_I_*L1*.tsv.gz | pigz -p 24 --fast --stdout > All_Intact_Samples_L1.tsv.gz
zcat *_I_*L2*.tsv.gz | pigz -p 12 --fast --stdout > All_Intact_Samples_L2.tsv.gz

mv All_Intact_Samples_L1.tsv.gz All_Intact_Samples_I_L1.tsv.gz
mv All_Intact_Samples_L2.tsv.gz All_Intact_Samples_I_L2.tsv.gz

######### BASH v
# "WK-1501-Test3_I_L1_FRAGSI_L1" "WK-1585-Citrate_I_L1_FRAGSI_L1" "WK-1585-Contrl_I_L1_FRAGSI_L1" "WK-1501-Test3_I_L2_FRAGSI_L2" "WK-1585-Citrate_I_L2_FRAGSI_L2" "WK-1585-Contrl_I_L2_FRAGSI_L2" 
exps=("All_Intact_Samples_I_L1" "All_Intact_Samples_I_L2")

==> run_macs3.sh <==
#!/bin/bash

# List of experiment names
exps=("CD28_L1_FRAGSCD28_L1" "CD28_L2_FRAGSCD28_L2" "All_Intact_Samples_I_L1" "All_Intact_Samples_I_L2")

# Export the macs3 command as a function so parallel can use it
run_macs3() {
  local EXP_name="$1"
  macs3 callpeak --treatment "${EXP_name}.tsv.gz" \
  --format BEDPE \
  --gsize mm \
  --nomodel \
  --llocal 50000 \
  --keep-dup all \
  --qvalue 0.01 \
  --call-summits \
  --outdir "${EXP_name}_macs3_Q01" \
  --name "${EXP_name}" \
  --verbose 2
}

# Export the function so GNU Parallel can use it
export -f run_macs3

# Run the commands in parallel using 20 jobs
parallel -j 24 run_macs3 ::: "${exps[@]}"

==> get_peak_numbers.sh <==
#!/bin/bash

# Output CSV file
output_file="summit_lines.csv"

# Write the header to the CSV file
echo "filename,number_of_lines" > "$output_file"

# Find all files in all subdirectories
find . -type f | while read -r file; do
  # Check if the file contains the word "summit"
  if grep -q "summit" "$file"; then
    # Count the number of lines in the file
    num_lines=$(wc -l < "$file")
    # Write the filename and line count to the CSV
    echo "$(basename "$file"),$num_lines" >> "$output_file"
  fi
done
######### BASH ^



import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the CSV file
csv_file = "summit_lines2.csv"
data = pd.read_csv(csv_file)

# Ensure the filename is treated as a string
data['filename'] = data['filename'].astype(str)

# Add columns to categorize files based on L1/L2 and content type (_I_ or CD28)
data['Type'] = data['filename'].apply(lambda x: 'Intact' if '_I_' in x else ('CD28' if 'CD28' in x else 'Other'))
data['Level'] = data['filename'].apply(lambda x: 'L1' if 'L1' in x else ('L2' if 'L2' in x else 'Unknown'))
# Extract everything before _I_ as a prefix if Type is 'I'
data['Label'] = data.apply(lambda row: row['filename'].split('_I_')[0]+"_Intact" if row['Type'] == 'Intact' else 'CD28', axis=1)

# Separate L1 and L2
data_L1 = data[data['Level'] == 'L1']
data_L2 = data[data['Level'] == 'L2']

# Set the aesthetic style for plots
sns.set(style="whitegrid")

# Function to create and save the plot
def create_and_save_plot(data_subset, level):
    # Plotting the data
    plt.figure(figsize=(10, 8))
    ax = sns.barplot(
        x='Label', 
        y='number_of_lines', 
        hue='Type', 
        data=data_subset, 
        palette='Set2'
    )
    
    # Add plot labels and title
    plt.title(f'{level} - Comparison of number of Peaks', fontsize=16)
    plt.xlabel(f'Time point', fontsize=14)
    plt.ylabel('Number of Peaks', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    
    # Save the plot to a PDF file
    output_file = f'{level}_peaks_comparison2.pdf'
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f'Saved plot: {output_file}')

# Create and save plots for L1 and L2
create_and_save_plot(data_L1, 'L1')
create_and_save_plot(data_L2, 'L2')



zcat All_Intact_Samples_I_L1.tsv.gz | wc -l
366506051
(scan) [ahrmad@electron Annot_spe_Peaks]$ zcat All_Intact_Samples_I_L2.tsv.gz | wc -l
1584466
(scan) [ahrmad@electron Annot_spe_Peaks]$ zcat CD28_L1_FRAGSCD28_L1.tsv.gz | wc -l
13316127
(scan) [ahrmad@electron Annot_spe_Peaks]$ zcat CD28_L2_FRAGSCD28_L2.tsv.gz | wc -l
4835380


import matplotlib.pyplot as plt
import snapatac2 as snap
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
from matplotlib_venn import venn3
import os
import csv

Experiment='Wouter21_SING_CB'
lum_ext = '_annot_Lum_newSigs.h5ad'

adata = sc.read_h5ad("Wouter21_SING_annot_All.h5ad")


adata_L1 = adata[np.isin(adata.obs['Annotation'],['L1'])].copy()

adata_L1_topTP = adata_L1[np.isin(adata_L1.obs['timePoint'],['Intact','Day28'])].copy()
