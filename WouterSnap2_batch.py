# Purpose:
#   Preprocess ATAC modality from 10x Multiome

#############################################
### Additional scripts to run before analysis vvv
#############################################

#############################################
###1 Peaks retrieval
#############################################
# to get BED file of peaks for each sample
# Used in the UNUSED Doublet Analysis based on peaks from 10x

##!/bin/bash
#
## Define the directory where you want to start the search
#root_directory="/mnt/ndata/daniele/wouter/Processed/CellRangerArc/"
#
#find "$root_directory" -type f -name "atac_peaks.bed" -path "*/outs/*" -exec sh -c '
#    folder_name=$(dirname "$0" | rev | cut -d/ -f2 | rev)
#    output_file="${folder_name}_atac_peaks.bed"
#    grep "^c" "$0" > "$output_file"
#' {} \;

#############################################
###2 Whitelists retrieval
#############################################
# to get the barcodes of cells as called by 10x pipeline
# Used in the import by SNAPATAC2

#!/bin/bash

## Set the directory where you have your subdirectories
#base_directory="/mnt/ndata/daniele/wouter/Processed/CellRangerArc"
#output_directory="/mnt/etemp/ahrmad/wouter/batch_ATAC"
#
## Loop over directories
#for directory in "$base_directory"/WK*; do
#    echo "$directory"
#    if [ -d "$directory" ]; then
#        csv_file=("$directory"/outs/per_barcode_metrics.csv)
#        summ_file=("$directory"/outs/summary.csv)
#        
#
#        mkdir -p "$output_directory/whitelists"
#        wl="$output_directory/whitelists/$(basename "$directory")"
#        wlLOG="$output_directory/whitelists/$(basename "$directory")_LOG"
#        awk -F ',' '$4 == 1 {print $1}' "$csv_file" > "$wl"
#
#        echo $(basename "$directory")",Total Barcodes,Cell Barcodes,ATAC Sequenced read pairs,GEX Sequenced read pairs" >> "${wlLOG}"
#        echo 0,$(wc -l "$csv_file" | cut -f1 -d" "),$(wc -l "$wl" | cut -f1 -d" "),$(cut -f23,45 -d',' "$summ_file"| sed -n '2 p') >> "${wlLOG}"
#
#        echo $(basename "$directory") whitelist done
#    fi
#done

#############################################
### Additional scripts to run before analysis ^^^
#############################################

from scipy.stats import median_abs_deviation
import snapatac2 as snap
import pandas as pd
import numpy as np
import os
print('snapatac2 version:',snap.__version__)

# Input files
path10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
pathRes = '/mnt/etemp/ahrmad/wouter/batch_ATAC'
exp10x = [d for d in os.listdir(path10x) if d.startswith('WK') and os.path.isdir(os.path.join(path10x, d))]

# QC filters
minFrags = 1000
minTSSe = 5

# Hard Threshold (cells are filtered out): maxFrags = max(MADs_maxFrags,Alt_maxFrags)
# Min Threshold (cells are labeled): maxFrags = Alt_maxFrags
Alt_maxFrags = 80000 
nmads = 6 

for exp in exp10x:
	print(f'{exp}')
	print(f'IMPORT AND FILTER')
	#Import ATAC fragments from 10x pipeline
	data = snap.pp.import_data(
		fragment_file=os.path.join(path10x,exp,'outs/atac_fragments.tsv.gz'),
		chrom_sizes=snap.genome.mm10,
		whitelist=os.path.join(pathRes,"whitelists",exp),
		sorted_by_barcode=False,min_num_fragments=0,
		tempdir=pathRes
		)

	#Filter cells
	snap.metrics.tsse(data, gene_anno=snap.genome.mm10, inplace=True)
	MADs_maxFrags = np.median(data.obs["n_fragment"]) + nmads * median_abs_deviation(data.obs["n_fragment"])
	maxFrags = max(MADs_maxFrags,Alt_maxFrags)
	snap.pp.filter_cells(data, min_counts=minFrags,max_counts=maxFrags, min_tsse=minTSSe, inplace=True)
	
	# QC metrics plots
	snap.metrics.frag_size_distr(data,add_key='frag_size_distr', inplace=True)
	snap.pl.tsse(data, min_fragment=0, width=750, height=600, interactive=False, show=False, out_file=exp+'_tsseFiltered.pdf')
	
	#Export fragments and anndata
	data.obs['Exp'] = pd.Categorical([exp]*data.n_obs)
	snap.ex.export_fragments(data, groupby='Exp', prefix='', suffix='.bed.gz')
	data.write(filename=f'{exp}_filt.h5ad')
	del data

# Run peak calling after 10x cellcalling and custom QC
# and IOM

#####################
############### BASH vvv
#####################
#!/bin/bash
# 1. Peak calling with MACS3
### PC.sh
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
        --outdir ${filename}_MACS_Q01 \
        --name $filename \
        --verbose 2
    echo "Peak calling done: $filename"
done
###

# 2. ATAC_IterativeOverlapMerging_v2.R
### PM.sh
#!/bin/bash
# Iterate through folders ending with "_MACS_Q01"
# and perform IterativeOverlapMerging of the peaks
# https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html
for folder in *_MACS_Q01; do
    if [ -d "$folder" ]; then
        # Get folder name without "_MACS_Q01"
        folder_name="${folder%_MACS_Q01}"
        echo "Peak Merging: $folder_name"

        # Change directory to the current folder
        cp IOM.R "$folder"
        cd "$folder"

        # Launch the R script with the specified arguments
        Rscript IOM.R "$(pwd)" 500 "${folder_name}_summits.bed" "${folder_name}_ITMPeaks.bed"
        echo "Peak Merging done: $folder_name"

        # Move back to the parent directory
        rm IOM.R
        cd ..
    fi
done
###

#Activate mac3 env, deac, activate R env
./PC.sh && mmde && mmac R && ./PM.sh

#####################
############### BASH ^^^
#####################

#Doublet Analysis based on tiles
	# for exp in exp10x:
		# print(f'DOUBLET ANALYSIS')
		# FiltData = snap.read(exp +'_filt.h5ad').to_memory()	
		# snap.pp.add_tile_matrix(FiltData, bin_size=1000, inplace=True)
		# snap.pp.select_features(FiltData, n_features=250000, inplace=True, blacklist="mm10-blacklist.v2.bed.gz") #
		# snap.pp.scrublet(FiltData, features='selected', n_comps=15, sim_doublet_ratio=2.0, expected_doublet_rate=0.1)	
		# doub = snap.pp.filter_doublets(FiltData, probability_threshold=0.7, inplace=False)
		# FiltData.obs['doublet_class'] = np.where(doub, 'singlet', 'doublet')
		# FiltData.write(filename=f'{exp}_filt_BINS.h5ad')

		# df = pd.DataFrame({'sample': exp, 'obs_names': FiltData.obs_names, 'doublet_score': FiltData.obs['doublet_probability'], 'doublet_class': FiltData.obs['doublet_class']})
		# df.to_csv(f'{exp}_doublet_scores_ATAC.csv', index=False)
		# del FiltData
		# print(f'{exp} DONE')

#Doublet Analysis based on peaks from 10x
	# for exp in exp10x:
		# print(f'DOUBLET ANALYSIS PEAKS')
		# pm = snap.read(exp +'_filt.h5ad').to_memory()
		# data = snap.pp.make_peak_matrix(pm, peak_file=f'{exp}_atac_peaks.bed')
		# data.obsm = pm.obsm.copy()
		# data.uns['reference_sequences'] = pm.uns['reference_sequences'].copy()
		# del pm
		# snap.pp.select_features(data, n_features=80000, inplace=True, blacklist="mm10-blacklist.v2.bed.gz") #
		# snap.pp.scrublet(data, features='selected', n_comps=15, sim_doublet_ratio=2.0, expected_doublet_rate=0.1)	
		# doub = snap.pp.filter_doublets(data, probability_threshold=0.7, inplace=False)
		# data.obs['doublet_class'] = np.where(doub, 'singlet', 'doublet')
		# data.write(filename=f'{exp}_filt_PEAKS.h5ad')

		# df = pd.DataFrame({'sample': exp, 'obs_names': data.obs_names, 'doublet_score': data.obs['doublet_probability'], 'doublet_class': data.obs['doublet_class']})
		# df.to_csv(f'{exp}_doublet_scores_ATACPeaks.csv', index=False)

#Doublet Analysis based on peaks from Macs3 (after 10x cellcalling and custom QC)
for exp in exp10x:
	print(f'DOUBLET ANALYSIS PEAKS SELF')

	# Read in QC'd Anndata file
	pm = snap.read(exp +'_filt.h5ad').to_memory()
	# Create cell-by-peak matrix from it with IOM'd Macs3 peak file
	data = snap.pp.make_peak_matrix(pm, peak_file=f'{exp}_MACS_Q01/{exp}_ITMPeaks.bed')
	# Copy obsm and ref seqs
	data.obsm = pm.obsm.copy()
	data.uns['reference_sequences'] = pm.uns['reference_sequences'].copy()
	del pm
	
	# Feature selection 
	snap.pp.select_features(data, n_features=80000, inplace=True, blacklist="mm10-blacklist.v2.bed.gz")
	# Doublet Detection
	snap.pp.scrublet(data, features='selected', n_comps=15, sim_doublet_ratio=2.0, expected_doublet_rate=0.1)	
	# Doublet Annotation
	doub = snap.pp.filter_doublets(data, probability_threshold=0.7, inplace=False)
	data.obs['doublet_class'] = np.where(doub, 'singlet', 'doublet')
	
	# Export Anndata
	data.write(filename=f'{exp}_filt_PEAKS_Self.h5ad')
	# Export CSV containing Doublet scores and annotation
	df = pd.DataFrame({'sample': exp, 'obs_names': data.obs_names, 'doublet_score': data.obs['doublet_probability'], 'doublet_class': data.obs['doublet_class']})
	df.to_csv(f'{exp}_doublet_scores_ATACPeaks_Self.csv', index=False)

	print(f'{exp} QC DONE and files created')

#######################################################################
#######################################################################
## Here normal analysis switches to a higher level (by merging samples)
#######################################################################
#######################################################################


print(f'Loading Doublet Information...')
# Doublet Analysis from RNA modality is required
# Create Whitelist Barcode Dictionary containing singlet cells by exp
BC_dict = {}
for exp in exp10x:
    atacdf = pd.read_csv(f'./csv/{exp}_doublet_scores_ATACPeaks_Self.csv')
    gexdf = pd.read_csv(f'./csv/{exp}_doublet_scores_GEX.csv')
    merged_df_all = atacdf.merge(gexdf, on='obs_names', how='inner')
    merged_df_all['doublet_class'] = 'WHATAMI'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet GEX Only'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Singlet ATAC Only'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
    
    # Retain QC passing cells (present in the CSV) 
    # that were called singlets by at least 1 modality
    merged_df_all_singlet = merged_df_all[merged_df_all['doublet_class'].str.contains('Singlet')]
    BC_dict[exp] = list(merged_df_all_singlet['obs_names'])

for exp in exp10x:

	# Import anndata from peak calling
	data = snap.read(f'{exp}_filt_PEAKS_Self.h5ad').to_memory()
	
	# Filter cells that are not in the Whitelist Barcode Dictionary	
	data = data[data.obs.index.isin(BC_dict[exp])].copy()

	# Analysis
	#Perform dimension reduction using the spectrum of the normalized graph Laplacian defined by pairwise similarity between cells
	snap.tl.spectral(data, weighted_by_sd=True,chunk_size=80000, features='selected', distance_metric='cosine', inplace=True)
	snap.tl.umap(data, use_rep='X_spectral', key_added='umap', random_state=0)
	#neighborhood graph of observations stored in adata using the method specified by method. The distance metric used is Euclidean.
	snap.pp.knn(data, n_neighbors=50, use_rep='X_spectral', method='kdtree')
	#Cluster cells using the Leiden algorithm [Traag18]
	snap.tl.leiden(data, resolution=0.5, objective_function='modularity', min_cluster_size=10)
	snap.pl.umap(data, color='leiden', height=500,interactive=False, show=False, out_file=exp+'_LeidenUMAP.pdf')

	# Export complete Anndata for the exp
	data.write(f'{exp}_postPeaks_Self.h5ad')

	# Create Gene matrix Anndata
	# Build Gene matrix
	gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.mm10)
	# Do some basic filtering 
	import scanpy as sc
	sc.pp.filter_genes(gene_matrix, min_cells= 5)
	sc.pp.normalize_total(gene_matrix)
	sc.pp.log1p(gene_matrix)
	sc.external.pp.magic(gene_matrix, solver="approximate")
	gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]
	# Export Gene matrix to Anndata
	gene_matrix.write(f'{exp}_gene_matPeaks_Self.h5ad')
	
	print(f'{exp} DONE')

# Annotation
import scanpy as sc
import matplotlib.pyplot as plt

snap_ext = '_gene_matPeaks_Self.h5ad'
input_h5ad = [d.removesuffix(snap_ext) for d in os.listdir('.') if d.endswith(snap_ext)]

for sample in input_h5ad:
	gene_matrix = sc.read_h5ad(sample + snap_ext)

	markerGenes = [
	    "Pbsn","Dpp4","Prom1",  #L1
	    "Ly6a","Tacstd2","Psca","Krt4", #L2
	    "Foxi1", #L3
	    "Krt8","Krt18" #All Luminal cells
	    ]
	
	markerGenes_L1 = ["Pbsn","Dpp4","Prom1"] #L1
	markerGenes_L2 = ["Ly6a","Tacstd2","Psca","Krt4"] #L2
	markerGenes_L3 = ["Foxi1"] #L3
	markerGenes_L = ["Krt8","Krt18"] # All Luminal

	sc.pl.umap(gene_matrix, use_raw=False, color=markerGenes_L1,save=exp+'_UMAP_GeneMatrix_L1.pdf',title=[exp+'_'+mark for mark in markerGenes_L1],show=False)
	sc.pl.umap(gene_matrix, use_raw=False, color=markerGenes_L2,save=exp+'_UMAP_GeneMatrix_L2.pdf',title=[exp+'_'+mark for mark in markerGenes_L2],show=False)
	sc.pl.umap(gene_matrix, use_raw=False, color=markerGenes_L3,save=exp+'_UMAP_GeneMatrix_L3.pdf',title=[exp+'_'+mark for mark in markerGenes_L3],show=False)
	sc.pl.umap(gene_matrix, use_raw=False, color=markerGenes_L,save=exp+'_UMAP_GeneMatrix_L.pdf',title=[exp+'_'+mark for mark in markerGenes_L],show=False)
	sc.pl.umap(gene_matrix, use_raw=False, color=markerGenes,save=exp+'_UMAP_GeneMatrix_ALL.pdf',title=[exp+'_'+mark for mark in markerGenes],show=False)


markerGenes = {
    "L1": ["Pbsn","Nkx3.1","CD26","Dpp4","CD59a","CD133","Prom1"],  #L1
    "L2": ["Sca1","Ly6a","Tacstd2","Trop2","Psca","Krt4","Claudin10"], #L2
    "L3": ["Foxi1","Atp6v1g3","Atp6b1b"], #L3
    "L": ["CD24a","Krt8","Krt18"] #All Luminal cells
    }
#SUPP
#luminal 1, Epcam, CD24a, Krt8, Krt18, Nkx3.1, Pbsn high
#luminal 2, Epcam, CD24a, Krt8, Krt18, Psca, Krt4, Tacst2, Ly6a
#luminal 3 (ionocyte), Epcam, CD24a, Krt8, Krt18, Foxi1, Atp6v1g3, Atp6b1b

for sample in input_h5ad:
	print(sample)
	adata = snap.read(sample + snap_ext).to_memory()
	adata.var.index = adata.var.index.str.upper()
	markerGenes_in_data = dict()

	for ct, markers in markerGenes.items():
	    markers_found = list()
	    for marker in markers:
	        if marker.upper() in adata.var.index:
	            markers_found.append(marker.upper())
	    markerGenes_in_data[ct] = markers_found
	for ct in list(markerGenes_in_data.keys()):
	    sc.pl.umap(
	        adata,
	        color=markerGenes_in_data[ct],
	        vmin=0,
	        vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
	        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
	        frameon=False,
	        cmap="plasma",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
	        return_fig=True
	    )
	    plt.suptitle(sample)
	    plt.savefig(sample+f'_UMAP_{ct.upper()}.pdf')
	    plt.close()

	sc.pl.dotplot(
	    adata,
	    groupby="leiden",
	    var_names=markerGenes_in_data,
	    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
	)
	plt.suptitle(sample)
	plt.subplots_adjust(bottom = 0.25)
	#plt.tight_layout()
	plt.savefig(sample+f'_DOTPLOT_LUMINAL.pdf')
	plt.close()















