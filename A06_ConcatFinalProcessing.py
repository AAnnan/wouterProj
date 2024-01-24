##Concatenate Samples

# Purpose:
#   Concatenate all samples in Experiment in ATAC
#   Perform the rest of the preprocessing

import snapatac2 as snap
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# To change according to samples
# SING for singulator, ENZ for enzymatic digestion
Experiments='Wouter21_ENZ'

# Input Files
DirATAC = '/mnt/etemp/ahrmad/wouter/batch_ATAC'
Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
qc_ext = '_filt.h5ad'
refDir = '/mnt/etemp/ahrmad/wouter/refs'
resDir = '/mnt/etemp/ahrmad/wouter/'
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

print(f'Concatenating using...')
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