# Purpose:
#   Merge Preprocessed samples in ATAC modality
#   Perform annotation

import snapatac2 as snap
import pandas as pd
import numpy as np
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
print('Anndata: ',ad.__version__,'Snap: ',snap.__version__)

# To change according to samples
# SING for singulator, ENZ for enzymatic digestion
Experiment='Wouter21_ENZ'

# Input Files
DirATAC = '/mnt/etemp/ahrmad/wouter/batch_ATAC'
Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
qc_ext = '_filt.h5ad'
refDir = '/mnt/etemp/ahrmad/wouter/refs'
All_Samples = [d for d in os.listdir(Dir10x) if d.startswith('WK')]

# Create a sample list based on the Experiment name
Samples = []
for sample in All_Samples:
	if '1350' in sample:
		if 'ENZ' in Experiment:
			Samples.append(sample)
	else:
		if 'SING' in Experiment:
			Samples.append(sample)

print(f'Loading Metadata...')
# Doublet Analysis from RNA modality is required
# Create Whitelist Barcode Dictionary containing singlet cells by exp
BC_dict = {}
for sample in Samples:
	atacdf = pd.read_csv(f'{refDir}/csv/{sample}_doublet_scores_ATACPeaks_Self.csv')
	gexdf = pd.read_csv(f'{refDir}/csv/{sample}_doublet_scores_GEX.csv')
	merged_df_all = atacdf.merge(gexdf, on='obs_names', how='inner')
	merged_df_all['doublet_class'] = 'WHATAMI'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet GEX Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Singlet ATAC Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
	
    # Retain QC passing cells (present in the CSV) 
    # that were called singlets by at least 1 modality
    merged_df_all_singlet = merged_df_all[merged_df_all['doublet_class'].str.contains('Singlet')]
	BC_dict[sample] = list(merged_df_all_singlet['obs_names'])

# Load metadata
tissueProv_dict = {}
for sample in Samples:
	if 'AP' in sample:
		prov='AP'
	elif 'LP' in sample:
		prov='LP'
	elif 'DP' in sample:
		prov='DP'
	elif 'VP' in sample:
		prov='VP'
	else:
		raise ValueError				
	tissueProv_dict[sample] = prov
isoMeth_dict = {}
for sample in Samples:
	if '1350' in sample:
		iso = 'Enzymatic Digestion'
	else:
		iso = 'Singulator'
	isoMeth_dict[sample] = iso
timePoint_dict = {}
for sample in Samples:
	if 'I' in sample:
		timeP='Intact'
	elif 'R3' in sample:
		timeP='RegenDay3'
	elif 'Regen_Day3' in sample:
		timeP='RegenDay3'
	elif 'RegenDay2' in sample:
		timeP='RegenDay2'
	elif 'RegenDay1' in sample:
		timeP='RegenDay1'
	elif 'Day28' in sample:
		timeP='Day28'
	else:
		raise ValueError
	timePoint_dict[sample] = timeP
mouseID_dict = {}
for sample in Samples:
	if 'I-1' in sample:
		mID = 'WK-1350-Intact-1'
	elif 'I-2' in sample:
		mID = 'WK-1350-Intact-2'
	elif 'R3-1' in sample:
		mID = 'WK-1350-R3-1'
	elif 'R3-2' in sample:
		mID = 'WK-1350-R3-2'
	else:
		mID = sample
	mouseID_dict[sample] = mID

# Concatenation
#Create concatenation list
adata_list = []
print(f'Loading AnnData...')
for sample in Samples:
	print(f'Reading {sample}')
	a = snap.read(os.path.join(DirATAC,sample+qc_ext)).to_memory()
	#a.obs['seqDate'] =  pd.Categorical([seqDate_dict[sample]]*a.n_obs)
	# Add metadata
	a.obs['Experiment'] = pd.Categorical([Experiment]*a.n_obs)
	a.obs['tissueProv'] =  pd.Categorical([tissueProv_dict[sample]]*a.n_obs)
	a.obs['timePoint'] =  pd.Categorical([timePoint_dict[sample]]*a.n_obs)
	a.obs['isoMeth'] =  pd.Categorical([isoMeth_dict[sample]]*a.n_obs)
	a.obs['mouseID'] =  pd.Categorical([mouseID_dict[sample]]*a.n_obs)

	# Filter cells that arent in the WL BC dict
	a = a[a.obs.index.isin(BC_dict[sample])].copy()
	# Rename the cells with their sample to avoid extrinsec duplicate
	a.obs.index = sample+ '_' + a.obs.index
	# Export
	a_fname = f'{sample}_qcTOREMOVE.h5ad'
	a.write(a_fname)
	# Append to list for concatenation
	adata_list.append((sample,a_fname))
	del a

print(f'Exporting fragments from all samples in {Experiment}...')
adata = snap.AnnDataSet(adata_list,filename=f'{Experiment}_ATAC.h5ad', add_key='sample')
adata.obs['Exp'] = [Experiment]*adata.n_obs
snap.ex.export_fragments(adata, groupby='Exp', prefix='', suffix='.bed.gz')
adata.close()
print(f'Launch "bash PC_IOM.sh {Experiment} <number of threads>"')

# Run peak calling after 10x cellcalling, custom QC and concat
# +IOM

#####################
############### BASH vvv
#####################

# 1. Peak calling with MACS3
### PC_IOM.sh
# (the .bed.gz is the fragment file exported in the previous step)
#!/bin/bash

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
#mm create -n PCIOM python=3.11
#mmac PCIOM
#pip install MACS3
#mm install r-base r-essentials r-data.table bioconductor-biocinstaller bioconductor-genomicranges
#OR
#mm env create --file PCIOM.yml

#Activate mac3 + R env before running
# bash PC_IOM.sh Wouter21_ENZ 20


#####################
############### BASH ^^^
#####################

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
print(f'Concatenating using...')
import scanpy as sc
import anndata as ad
print('Anndata: ',ad.__version__,'Scanpy: ',sc.__version__)

qc_ext='_qcTOREMOVE.h5ad'
#Concatenation
# building concat list
adata_list = []
for sample in Samples:
    print(f'Reading {sample}')
    a = sc.read_h5ad(sample+qc_ext)
    del a.obsm
    adata_list.append(a)
    del a
print(f'Concatenating...')
adata = ad.concat(adata_list, join='inner', merge='same',label='batch',keys=Samples,index_unique='_')
print(f'Writing {Experiment}_allQC.h5ad')
adata.write(Experiment + '_allQC.h5ad')

print(f'Analysis...')
data = snap.read(Experiment + '_allQC.h5ad').to_memory()

b = snap.pp.import_data(fragment_file=f'{Experiment}.bed.gz',
	chrom_sizes=snap.genome.mm10,
	sorted_by_barcode=False,min_num_fragments=0,
	tempdir='.')

# Reindexing
data.obs.index = [el.split('_WK')[0] for el in data.obs.index]
data.obs = data.obs.reindex(index=b.obs.index)
# Check if reindex was done properly
print('reindex done properly?')
print(np.all(data.obs.index == b.obs.index))
print(data.obs.index.equals(b.obs.index))
print(data.obs['batch'].value_counts(dropna=False))
print('------------')

# Get fragments,ref from b
data.obsm = b.obsm.copy()
data.uns['reference_sequences'] = b.uns['reference_sequences'].copy()

# Build cell-by-peak matrix and store in pm
pm = snap.pp.make_peak_matrix(data, peak_file=f'{Experiment}_ITMPeaks.bed')
# Copy fragments,ref from data
pm.obsm = data.obsm.copy()
pm.uns['reference_sequences'] = data.uns['reference_sequences'].copy()
del data,b

# Calculate FRiP 
snap.metrics.frip(pm, {"FRiP": f'{Experiment}_ITMPeaks.bed.gz'}, inplace=True)


# Feature Selection
#import urllib.request
#urllib.request.urlretrieve("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz", "mm10-blacklist.v2.bed.gz")
snap.pp.select_features(pm, n_features=120000, inplace=True, blacklist=f"{refDir}/mm10-blacklist.v2.bed.gz") #

#Perform dimension reduction using the spectrum of the normalized graph Laplacian defined by pairwise similarity between cells
snap.tl.spectral(pm, weighted_by_sd=True, chunk_size=80000, features='selected', distance_metric='cosine', inplace=True)
snap.tl.umap(pm, use_rep='X_spectral', key_added='umap', random_state=None)

#neighborhood graph of observations stored in data using the method specified by method. The distance metric used is Euclidean.
snap.pp.knn(pm, n_neighbors=50, use_rep='X_spectral', method='kdtree')

#Cluster cells using the Leiden algorithm [Traag18]
for resLeiden in [.25,.5,1,1.5,2]:
	print(f'Leiden clustering at {resLeiden} resolution')
	snap.tl.leiden(pm, resolution=resLeiden, key_added=f"leiden_res{resLeiden}")
	snap.pl.umap(pm, color=f"leiden_res{resLeiden}", height=500,interactive=False, show=False, out_file=f"{Experiment}_leiden_res{resLeiden}.pdf")
	sc.pl.umap(pm,color=f"leiden_res{resLeiden}",legend_loc="on data",save=f"{Experiment}_leiden_res{resLeiden}_sc.pdf",title=Experiment,show=False)

# Export Final peak matrix anndata
pm.write(f'{Experiment}_Post.h5ad')

#pm = snap.read(Experiment + '_Post.h5ad').to_memory()

# save OTHER UMAP PLOTS
sc.pl.umap(pm,color=["n_fragment", "frac_mito", "tsse"],save=Experiment+'_UMAP_QC',show=False)
sc.pl.umap(pm,color=["batch"],save=Experiment+'_batch',title=f'{Experiment}_ATAC',show=False)
sc.pl.umap(pm,color=["timePoint"],save=Experiment+'_timePoint',title=f'{Experiment}_ATAC',show=False)
sc.pl.umap(pm,color=["isoMeth"],save=Experiment+'_isoMeth',title=f'{Experiment}_ATAC',show=False)
sc.pl.umap(pm,color=["mouseID"],save=Experiment+'_mouseID',title=f'{Experiment}_ATAC',show=False)
#sc.pl.umap(pm,color=["seqDate"],save=Experiment+'_seqDate',title=f'{Experiment}_ATAC',show=False)
sc.pl.umap(pm,color=["tissueProv"],save=Experiment+'_tissueProv',title=f'{Experiment}_ATAC',show=False)

# Create Gene matrix Anndata
# Build Gene matrix
gene_matrix = snap.pp.make_gene_matrix(pm, snap.genome.mm10)
# Do some basic filtering 
sc.pp.filter_genes(gene_matrix, min_cells= 5)
sc.pp.normalize_total(gene_matrix)
sc.pp.log1p(gene_matrix)
sc.external.pp.magic(gene_matrix, solver="approximate")
gene_matrix.obsm["X_umap"] = pm.obsm["X_umap"]
# Export Gene matrix to Anndata
gene_matrix.write(f'{Experiment}_GeneMat.h5ad')

# Annotation
import scanpy as sc
import matplotlib.pyplot as plt

snap_ext = '_GeneMat.h5ad'

adata = sc.read_h5ad(Experiment + snap_ext)
#adata = gene_matrix

#UMAP GENE COLORS
#adata.var.index = adata.var.index.str.upper()
#adata.raw.var_names = adata.raw.var_names.str.upper()


#Ori p1
markerGenes = {
    "L1": ["Pbsn","Nkx3-1","CD26","Dpp4","CD59a","CD133","Prom1","DPP4"],  #L1
    "L2": ["Sca1","Ly6a","Tacstd2","Trop2","Psca","Krt4","Claudin10"], #L2
    "L3": ["Foxi1","Atp6v1g3","Atp6b1b"], #L3
    "L": ["CD24a","Krt8","Krt18","Pax2"] #All Luminal cells
    }
#SUPP p.12
#luminal 1, Epcam, CD24a, Krt8, Krt18, Nkx3.1, Pbsn high
#luminal 2, Epcam, CD24a, Krt8, Krt18, Psca, Krt4, Tacst2, Ly6a
#luminal 3 (ionocyte), Epcam, CD24a, Krt8, Krt18, Foxi1, Atp6v1g3, Atp6b1b

markerGenes_in_data = dict()
for ct, markers in markerGenes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    markerGenes_in_data[ct] = markers_found

for ct in list(markerGenes_in_data.keys()):
    print(f"{ct.upper()}")  # print cell subtype name
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
    plt.suptitle(Experiment)
    plt.savefig(Experiment+f'_UMAP_atac_{ct.upper()}.pdf')
    plt.close()

#DOTPLOT
sc.pl.dotplot(
    adata,
    groupby="leiden_res0.5",
    var_names=markerGenes_in_data,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
plt.suptitle(Experiment)
plt.subplots_adjust(bottom = 0.25)
plt.tight_layout()
plt.savefig(Experiment+f'_atac_DOTPLOT_LUMINAL.pdf')
plt.close()

#Celltyping based on Paper genes 
sc.tl.score_genes(adata, [m for m in markerGenes['L1']], ctrl_size=50, n_bins=25, score_name='L1_genes_Paper')
sc.tl.score_genes(adata, [m for m in markerGenes['L2']], ctrl_size=50, n_bins=25, score_name='L2_genes_Paper')
sc.tl.score_genes(adata, [m for m in markerGenes['L3']], ctrl_size=50, n_bins=25, score_name='L3_genes_Paper')
sc.tl.score_genes(adata, [m for m in markerGenes['L']], ctrl_size=50, n_bins=25, score_name='L_genes_Paper')

sc.pl.umap(adata,color=["L1_genes_Paper"],save=Experiment+'_atac_L1_genes_fromMain',title=f'{Experiment}_atac_L1 sig paper',show=False)
sc.pl.umap(adata,color=["L2_genes_Paper"],save=Experiment+'_atac_L2_genes_fromMain',title=f'{Experiment}_atac_L2 sig paper',show=False)
sc.pl.umap(adata,color=["L3_genes_Paper"],save=Experiment+'_atac_L3_genes_fromMain',title=f'{Experiment}_atac_L3 sig paper',show=False)
sc.pl.umap(adata,color=["L_genes_Paper"],save=Experiment+'_atac_L_genes_fromMain',title=f'{Experiment}_atac_L sig paper',show=False)

#Celltyping based on Supp gene lists
gene_set_table = f'{refDir}/table_s8_summary.txt'
gene_set_df = pd.read_csv(gene_set_table, sep='\t')

for gene_set in gene_set_df.columns:
	ctrl = 100
	print(f'{gene_set} S8')
	genes = gene_set_df[gene_set]
	sc.tl.score_genes(adata, genes, ctrl_size=ctrl, n_bins=25, score_name=gene_set)
	sc.pl.umap(adata,color=gene_set,save=f'{Experiment}_atac_{gene_set}_SuppS8_ctrl{ctrl}',title=f'{Experiment}_{gene_set}',show=False)

main_ct_marker_list = {"Epithelial": ["Epcam"], "Immune": ["Ptprc"]}  # present

sub_ct_marker_list = {
    "Bcells": ["Cd19", "Ms4a1"],
    "Tcells": ["Cd3g", "Cd3d", "Cd3e", "Cd247"],
    "macrophages": ["Cd14", "Aif1"],
    "endothelium": ["Pecam1", "Vwf"],
    "lymphatic_endothelium": ["Pecam1", "Prox1"],
    "glia": ["Sox10"],
    "myofibroblast": ["Acta2", "Myh11", "Rspo3"],
    "smooth_muscle": ["Acta2", "Notch3"],
    "mesenchymal_1": ["Col5a2", "Lama2", "Zeb1", "Wnt2", "Wnt6", "Wnt10a", "Rorb"],
    "mesenchymal_2": ["Col5a2", "Lama2", "Zeb1", "Sult1e1", "Fgf10", "Rspo1"],
    "basal": ["Epcam", "Krt5", "Krt14", "Trp63"],
    "seminal_vesicle_basal": ["Epcam", "Pax2", "Krt5", "Krt14", "Trp63", "Calml3"],
    "luminal_1": ["Epcam", "Cd24a", "Krt8", "Krt18", "Nkx3-1", "Pbsn", "PROM1", "CD26/DPP4"],
    "luminal_2": ["Epcam", "Cd24a", "Krt8", "Krt18", "Psca", "Krt4", "Tacstd2", "Ly6a"],
    "luminal_3": ["Epcam", "Cd24a", "Krt8", "Krt18", "Foxi1", "Atp6v1g3", "Atp6v1b1"],
    "seminal_vesicle_luminal": ["Epcam", "Pax2", "Krt8", "Krt18", "Pate4"],
    "seminal_vesicle_ionocyte": ["Epcam", "Pax2", "Foxi1"],
}  # already alias-corrected

# big absences from 'measured_genes_all': none
# big absences from 'measured_genes' (Experiment-specific): Cd19, Cd3*, Cd31/Pecam1, Sox10, some Wnt, Cd24a, ("Foxi1", "Atp6v1g3", "Atp6v1b1") discriminating the luminal_3

for gene_set in sub_ct_marker_list.keys():
	print(f'{gene_set} DanDict')
	genes = pd.Series(sub_ct_marker_list[gene_set])
	sc.tl.score_genes(adata, genes, ctrl_size=50, n_bins=25, score_name=gene_set)
	sc.pl.umap(adata,color=gene_set,save=f'{Experiment}_atac_{gene_set}_sub_ct_marker',title=f'{Experiment}_{gene_set}_sub_ct_marker',show=False)

for gene_set in main_ct_marker_list.keys():
	print(f'{gene_set} DanDict')
	genes = pd.Series(main_ct_marker_list[gene_set])
	sc.tl.score_genes(adata, genes, ctrl_size=50, n_bins=25, score_name=gene_set)
	sc.pl.umap(adata,color=gene_set,save=f'{Experiment}_atac_{gene_set}_main_ct_marker',title=f'{Experiment}_{gene_set}_main_ct_marker',show=False)






