from scipy.sparse import csr_matrix, issparse
import scanpy as sc
import pandas as pd
import numpy as np
import os


#import anndata
#import matplotlib.pyplot as plt
#import scvi
#import torch
#from scvi.model.utils import mde

#sc.set_figure_params(figsize=(5, 5), frameon=False)
#torch.set_float32_matmul_precision("high")
#refDir = '/scratch/aannan'

refDir = '/mnt/etemp/ahrmad/wouter/refs'

# File Paths
k20_log2TP10k_path = f'{refDir}/Karthaus2020_data/GSE146811_mmProstate10x_timecourse_log2TP10k.h5'
k20_rawCounts_path = f'{refDir}/Karthaus2020_data/GSE146811_mmProstate10x_timecourse_rawCount.h5'
k20_full_sample_final_path = f'{refDir}/Karthaus2020_data/GSE146811_mmProstate10x_full_sample_final.tsv'

#w21_SING_rna_path = f'{refDir}/Karthaus2020_data/Wouter21_SING_CB_postToMap.h5ad'
#w21_ENZ_rna_path = f'{refDir}/Karthaus2020_data/Wouter21_ENZ_CB_postToMap.h5ad'

w21_SING_rna_path = f'/mnt/etemp/ahrmad/wouter/ENZ_RNA_CB/Wouter21_ENZ_CB_post.h5ad'
w21_ENZ_rna_path = f'/mnt/etemp/ahrmad/wouter/SING_RNA_CB/Wouter21_SING_CB_post.h5ad'

# Import and prepare Adata
k20_log2TP10k = sc.read_10x_h5(k20_log2TP10k_path)
k20_log2TP10k.var['gene_names'] = k20_log2TP10k.var.index
k20_log2TP10k.var.index = k20_log2TP10k.var['gene_ids']
k20_log2TP10k.layers["counts"] = k20_log2TP10k.X

k20_rawCounts = sc.read_10x_h5(k20_rawCounts_path)
k20_rawCounts.var['gene_names'] = k20_rawCounts.var.index
k20_rawCounts.var.index = k20_rawCounts.var['gene_ids']
k20_rawCounts.layers["counts"] = k20_rawCounts.X

k20_full_sample_final = pd.read_csv(k20_full_sample_final_path,sep='\t',header=0,index_col=0)

#for col in k20_full_sample_final.columns[2:]:
#	print(col)
#	print(np.unique(k20_full_sample_final[col], return_counts=True))
#	print('\n\n\n')

#k20_full_sample_final['predTypeInt'].unique().shape

w21_SING_rna = sc.read_h5ad(w21_SING_rna_path)
w21_SING_rna.layers["counts"] = w21_SING_rna.X
w21_SING_rna.var['gene_names'] = w21_SING_rna.var.index
w21_SING_rna.var.index = w21_SING_rna.var['gene_ids']

w21_ENZ_rna = sc.read_h5ad(w21_ENZ_rna_path)
w21_ENZ_rna.layers["counts"] = w21_ENZ_rna.X
w21_ENZ_rna.var['gene_names'] = w21_ENZ_rna.var.index
w21_ENZ_rna.var.index = w21_ENZ_rna.var['gene_ids']

# Add missing genes with counts 0 to map
missing_genes_enz = [
    gene_id
    for gene_id in k20_log2TP10k.var.index
    if gene_id not in w21_ENZ_rna.var.index
]

missing_genes_sing = [
    gene_id
    for gene_id in k20_log2TP10k.var.index
    if gene_id not in w21_SING_rna.var.index
]

missing_gene_w21_ENZ_rna = sc.AnnData(
    X=csr_matrix(np.zeros(shape=(w21_ENZ_rna.n_obs, len(missing_genes_enz))), dtype="float32"),
    obs=w21_ENZ_rna.obs.iloc[:, :1],
    var=k20_log2TP10k.var.loc[missing_genes_enz, :])
missing_gene_w21_ENZ_rna.layers["counts"] = missing_gene_w21_ENZ_rna.X

missing_gene_w21_SING_rna = sc.AnnData(
    X=csr_matrix(np.zeros(shape=(w21_SING_rna.n_obs, len(missing_genes_sing))), dtype="float32"),
    obs=w21_SING_rna.obs.iloc[:, :1],
    var=k20_log2TP10k.var.loc[missing_genes_sing, :])
missing_gene_w21_SING_rna.layers["counts"] = missing_gene_w21_SING_rna.X


w21_ENZ_rna_augmented = sc.concat(
    [w21_ENZ_rna, missing_gene_w21_ENZ_rna],
    axis=1,
    join="outer",
    index_unique=None,
    merge="unique",
)
w21_SING_rna_augmented = sc.concat(
    [w21_SING_rna, missing_gene_w21_SING_rna],
    axis=1,
    join="outer",
    index_unique=None,
    merge="unique",
)

w21_ENZ_rna_augmented = w21_ENZ_rna_augmented[
    :, k20_log2TP10k.var.index
].copy()

w21_SING_rna_augmented = w21_SING_rna_augmented[
    :, k20_log2TP10k.var.index
].copy()


assert (k20_rawCounts.var.index == k20_log2TP10k.var.index).all()
assert (w21_ENZ_rna_augmented.var.index == k20_log2TP10k.var.index).all()
assert (w21_SING_rna_augmented.var.index == k20_log2TP10k.var.index).all()

# Change gene ids back to gene names
#w21_ENZ_rna_augmented.var["gene_ids"] = w21_ENZ_rna_augmented.var.index
#w21_ENZ_rna_augmented.var.set_index("gene_names", inplace=True)
#
#w21_SING_rna_augmented.var["gene_ids"] = w21_SING_rna_augmented.var.index
#w21_SING_rna_augmented.var.set_index("gene_names", inplace=True)


w21_ENZ_rna_augmented.obs['samples'] = w21_ENZ_rna_augmented.obs.batch 
w21_ENZ_rna_augmented.obs.batch = 'ENZ'
#w21_ENZ_rna_augmented.obs.batch.unique()
w21_SING_rna_augmented.obs['samples'] = w21_SING_rna_augmented.obs.batch 
w21_SING_rna_augmented.obs.batch = 'SING'
#w21_SING_rna_augmented.obs.batch.unique()




#############
# INGEST
w21_ENZ_rna = w21_ENZ_rna[:,w21_ENZ_rna.var.index.isin(k20_log2TP10k.var.index)]
w21_SING_rna = w21_SING_rna[:,w21_SING_rna.var.index.isin(k20_log2TP10k.var.index)]
k20_log2TP10k = k20_log2TP10k[:,k20_log2TP10k.var.index.isin(w21_ENZ_rna.var.index)]

w21_ENZ_rna.var = w21_ENZ_rna.var.reindex(index=sorted(w21_ENZ_rna.var.index))
w21_SING_rna.var = w21_SING_rna.var.reindex(index=sorted(w21_SING_rna.var.index))
k20_log2TP10k.var = k20_log2TP10k.var.reindex(index=sorted(k20_log2TP10k.var.index))


assert (w21_ENZ_rna.var.index == k20_log2TP10k.var.index).all()
assert (w21_SING_rna.var.index == k20_log2TP10k.var.index).all()

adata_ref = k20_log2TP10k.copy()

sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)

assert (adata_ref.obs.index == k20_full_sample_final.index).all()
for col in k20_full_sample_final.columns[2:]:
	adata_ref.obs[col] = k20_full_sample_final[col]

sc.tl.ingest(w21_ENZ_rna, adata_ref, obs='predTypeInt')
sc.tl.ingest(w21_SING_rna, adata_ref, obs='predTypeInt')

w21_ENZ_rna.obs['predTypeInt']


sc.set_figure_params(figsize=(5, 5))
sc.settings.set_figure_params(
    dpi_save=200,
    fontsize=12,
    facecolor="white",
    frameon=False,
)
sc.pl.umap(w21_ENZ_rna,color='predTypeInt',save=f'Wouter21_ENZ_LABELED.png',title=f'ENZ Label Transfer Karthaus 2020',show=False)
sc.pl.umap(w21_SING_rna,color='predTypeInt',save=f'Wouter21_SING_LABELED.png',title=f'SING Label Transfer Karthaus 2020',show=False)


w21_SING_rna_path = f'/mnt/etemp/ahrmad/wouter/ENZ_RNA_CB/Wouter21_ENZ_CB_post.h5ad'
w21_ENZ_rna_path = f'/mnt/etemp/ahrmad/wouter/SING_RNA_CB/Wouter21_SING_CB_post.h5ad'

tt = sc.read_h5ad(w21_SING_rna_path)
sc.pl.umap(tt,color='timePoint',save=f'TEST.png',title=f'TT',show=False)



############
##scanVI

# Prepare scanVI
scvi.settings.seed = 0
print("scvi-tools version:", scvi.__version__)

assert (k20_rawCounts.obs.index == k20_full_sample_final.index).all()
for col in k20_full_sample_final.columns[2:]:
	k20_rawCounts.obs[col] = k20_full_sample_final[col]


#np.unique(k20_log2TP10k.obs['predType'], return_counts=True)

# Train model
k20_rawCounts.layers["counts"] = k20_rawCounts.X.copy()
sc.pp.normalize_total(k20_rawCounts, target_sum=1e4)
sc.pp.log1p(k20_rawCounts)
k20_rawCounts.raw = k20_rawCounts  # keep full dimension safe

#sc.pp.highly_variable_genes(k20_rawCounts,flavor="seurat_v3",n_top_genes=4000,layer="counts",subset=True)
scvi.model.SCVI.setup_anndata(k20_rawCounts, layer="counts")
scvi_model = scvi.model.SCVI(k20_rawCounts, n_layers=2, n_latent=30)
scvi_model.train()

SCVI_LATENT_KEY = "X_scVI"
k20_rawCounts.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()


#Transfer
SCANVI_CELLTYPE_KEY = "celltype_scanvi"
#ENZ

w21_ENZ_rna_augmented.obs[SCANVI_CELLTYPE_KEY] = "Unknown"

scvi.model.SCVI.setup_anndata(w21_ENZ_rna_augmented, layer="counts")

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=w21_ENZ_rna_augmented,
    unlabeled_category="Unknown",
    labels_key=SCANVI_CELLTYPE_KEY,
)

scanvi_model.train(max_epochs=30, n_samples_per_label=100)

SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTION_KEY = "C_scANVI"

w21_ENZ_rna_augmented.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(w21_ENZ_rna_augmented)
w21_ENZ_rna_augmented.obs[SCANVI_PREDICTION_KEY] = scanvi_model.predict(w21_ENZ_rna_augmented)

SCANVI_MDE_KEY = "X_mde_scanvi"
w21_ENZ_rna_augmented.obsm[SCANVI_MDE_KEY] = mde(w21_ENZ_rna_augmented.obsm[SCANVI_LATENT_KEY])
w21_ENZ_rna_augmented.write('Wouter21_ENZ_CB_scVI.h5ad')

#SING
w21_SING_rna_augmented.obs[SCANVI_CELLTYPE_KEY] = "Unknown"

scvi.model.SCVI.setup_anndata(w21_SING_rna_augmented, layer="counts")

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=w21_SING_rna_augmented,
    unlabeled_category="Unknown",
    labels_key=SCANVI_CELLTYPE_KEY,
)

scanvi_model.train(max_epochs=30, n_samples_per_label=100)

SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTION_KEY = "C_scANVI"

w21_SING_rna_augmented.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(w21_SING_rna_augmented)
w21_SING_rna_augmented.obs[SCANVI_PREDICTION_KEY] = scanvi_model.predict(w21_SING_rna_augmented)

SCANVI_MDE_KEY = "X_mde_scanvi"
w21_SING_rna_augmented.obsm[SCANVI_MDE_KEY] = mde(w21_SING_rna_augmented.obsm[SCANVI_LATENT_KEY])
w21_SING_rna_augmented.write('Wouter21_SING_CB_scVI.h5ad')











#!/bin/bash

##############################
## SCRIPT RUNNING ON UNIL GPUs
##############################
#! /bin/bash

## Resource Allocation
#SBATCH --time=1-00:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32

## job metadata
#SBATCH --job-name="Label Transfer"
#SBATCH --mail-user=ahrmad.annan@unil.ch
#SBATCH --mail-type=end,fail
set -e

source /scratch/aannan/micromamba/etc/profile.d/micromamba.sh

micromamba activate label
srun python label.py
micromamba deactivate 

#sbatch go.sh

mm install python==3.11 pip ipython
mm install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
pip install scvi-tools scanpy scikit-misc scipy pandas numpy anndata matplotlib
pip install jax https://storage.googleapis.com/jax-releases/cuda12/jaxlib-0.4.23+cuda12.cudnn89-cp311-cp311-manylinux2014_x86_64.whl


pip uninstall jax jax-cuda12-pjrt jax-cuda12-plugin jaxlib



