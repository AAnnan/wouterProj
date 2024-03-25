import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix, issparse
import anndata as ad
import os

import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi_save=200,
    fontsize=12,
    facecolor="white",
    frameon=False,
)
print('Anndata: ',ad.__version__,'Scanpy: ',sc.__version__)


### TO CHECK WHEN CHANGING SAMPLES ###
Experiment='Wouter21_ENZ_CB'
### TO CHECK WHEN CHANGING SAMPLES ###
DirRNA = '/mnt/etemp/ahrmad/wouter/RNA_CB'
refDir = '/mnt/etemp/ahrmad/wouter/refs'
if not os.path.exists(refDir):
    raise ValueError

#Naming
qc_ext = '_qc.h5ad'
post_ext = '_post.h5ad'
atac_doub_csv_ext = '_doublet_scores_ATACPeaks_Self.csv'
gex_doub_csv_ext = '_doublet_scores_CB_GEX.csv'


Samps = [s.split(qc_ext)[0] for s in os.listdir(DirRNA) if s.startswith('WK')]
Samples = []
if 'ENZ' in Experiment:
    Samples.extend([s for s in Samps if '1350' in s])
elif 'SING' in Experiment:
    Samples.extend([s for s in Samps if '1350' not in s])
else:
    raise ValueError

print(f'Loading Samples Metadata...')

seqDate_dict = {}
with open(f'{refDir}/sequencing_batches.txt', 'r') as seqbatchfile:
    next(seqbatchfile)
    for line in seqbatchfile:
        key, value = line.strip().split('\t')
        seqDate_dict[key] = value

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

isoMeth_dict = {}
for sample in Samples:
    if '1350' in sample:
        iso = 'Enzymatic Digestion'
    else:
        iso = 'Singulator'
    isoMeth_dict[sample] = iso

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


### Concatenation
adata_list = []
print(f'Loading AnnData...')
for sample in Samples:
    print(f'\tReading {sample}')
    sample_adata = sc.read_h5ad(os.path.join(DirRNA,sample+qc_ext))
    sample_adata.obs['tissueProv'] = tissueProv_dict[sample]
    sample_adata.obs['seqDate'] = seqDate_dict[sample]
    sample_adata.obs['timePoint'] = timePoint_dict[sample]
    sample_adata.obs['isoMeth'] = isoMeth_dict[sample]
    sample_adata.obs['mouseID'] = mouseID_dict[sample]
    
    #Filter doublets per sample
    atacdf = pd.read_csv(f'{refDir}/csv/{sample}{atac_doub_csv_ext}')
    gexdf = pd.read_csv(f'{refDir}/csv/{sample}{gex_doub_csv_ext}')
    print(sample)
    print(f'Pre-merging: ATAC: {atacdf.shape[0]} GEX: {gexdf.shape[0]} cells')
    mergDF = atacdf.merge(gexdf, on='obs_names', how='inner')
    mergDF['doublet_class'] = 'WHATAMI'
    mergDF.loc[(mergDF['doublet_class_x'] == 'singlet') & (mergDF['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
    mergDF.loc[(mergDF['doublet_class_x'] == 'doublet') & (mergDF['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet GEX Only'
    mergDF.loc[(mergDF['doublet_class_x'] == 'singlet') & (mergDF['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Singlet ATAC Only'
    mergDF.loc[(mergDF['doublet_class_x'] == 'doublet') & (mergDF['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
    
    # Check if the doublet_class column contains 'WHATAMI'
    if any(mergDF['doublet_class'].str.contains('WHATAMI')):
        raise ValueError(f"Some cells have doublet status unassigned")
    print(f'Post-merging: {mergDF.shape[0]} cells')
    
    #Select only singlets
    mergDF_singlet = mergDF[mergDF['doublet_class'] == 'Singlet Only']

    #mergDF_singlet = mergDF.copy() #-> for doublet comparison between AMU and SCRUBLET

    print(f'Singlets: {mergDF_singlet.shape[0]} cells')

    sample_adata = sample_adata[np.isin(sample_adata.obs.index, mergDF_singlet['obs_names'])].copy()
    
    mergDF_singlet.index = mergDF_singlet['obs_names']
    mergDF_singlet_sorted = mergDF_singlet.reindex(sample_adata.obs.index)
    assert list(sample_adata.obs_names) == list(mergDF_singlet_sorted['obs_names']),'Adata and dblt CSV not in same order'
    sample_adata.obs['doublet_class'] = mergDF_singlet_sorted['doublet_class']

    adata_list.append(sample_adata)


print(f'Concatenating...')
adata = ad.concat(adata_list, join='inner', merge='same',label='batch',keys=Samples,index_unique='_')

### SCRAN Normalization
print('Starting SCRAN Normalization\nlog1p with Scran estimated size factors')
#log1p with Scran estimated size factors
### Feature selection and Dimensionality Reduction
%R library(scran)
%R library(BiocParallel)
%R library(scry)
xTopDeviantGenes = 4000

# Preliminary clustering for differentiated normalisation
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups")
data_mat = adata_pp.X.T

# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()

ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata_pp.obs["groups"]
del adata_pp
%R size_factors = sizeFactors(computeSumFactors(SingleCellExperiment(list(counts=data_mat)),clusters = input_groups,min.mean = 0.1,BPPARAM = MulticoreParam()))

adata.obs["size_factors"] = %Rget size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("log1p with Scran estimated size factors")
fig.suptitle(Experiment)
if not os.path.exists('figures'):
    os.makedirs('figures')
plt.savefig(os.getcwd()+'/figures/'+sample+'_Norm.pdf')
plt.close()

print('SCRAN Normalization DONE')

### Feature selection and Dimensionality Reduction
print('Starting Feature selection and Dimensionality Reduction')
print('Deviance')
ro.globalenv["adata"] = adata

# Feature selection 
%R sce = devianceFeatureSelection(adata, assay="X")

binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T

idx = binomial_deviance.argsort()[-xTopDeviantGenes:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance

#compute the mean and dispersion for each gene accross all cells
sc.pp.highly_variable_genes(adata, layer="scran_normalization")

sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5)
plt.title(Experiment)
plt.savefig(os.getcwd()+'/figures/'+sample+'_Disp.pdf')
plt.close()

print(f'{Experiment} Feature Selection with Deviance DONE')

adata.layers["No_normalization"] = adata.X.copy()
adata.X = adata.layers["scran_normalization"]

#PCA
# setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_variance_ratio(adata, log=True,save=Experiment+'_PCA_Variance',show=False)

sc.pp.neighbors(adata, use_rep="X_pca",n_neighbors=30, n_pcs=20)
sc.tl.umap(adata)

for resLeiden in [.25,.5,1,1.25,1.5,2,2.5,3]:
    print(f'Leiden clustering at {resLeiden} resolution')
    sc.tl.leiden(adata, key_added=f"leiden_res{resLeiden}", resolution=resLeiden)
    sc.pl.umap(adata,color=f"leiden_res{resLeiden}",legend_loc="on data",save=f"{Experiment}_leiden_res{resLeiden}",title=Experiment,show=False)

#OTHER UMAP PLOTS
#QC
sc.pl.umap(adata,color=["scDblFinder_score", "scDblFinder_class"],save=Experiment+'_UMAP_Dblt_QC',show=False)
sc.pl.umap(adata,color=["doublet_class"],save=Experiment+'_UMAP_Dblt_QC',show=False)
sc.pl.umap(adata,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'_UMAP_Counts_QC',show=False)
sc.pl.umap(adata,color=["pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],save=Experiment+'_UMAP_PctCounts_QC',show=False)
sc.pl.umap(adata,color=["batch"],save=Experiment+'_batch',title=Experiment,show=False)
sc.pl.umap(adata,color=["timePoint"],save=Experiment+'_timePoint',title=Experiment,show=False)
sc.pl.umap(adata,color=["isoMeth"],save=Experiment+'_isoMeth',title=Experiment,show=False)
sc.pl.umap(adata,color=["mouseID"],save=Experiment+'_mouseID',title=Experiment,show=False)
sc.pl.umap(adata,color=["seqDate"],save=Experiment+'_seqDate',title=Experiment,show=False)
sc.pl.umap(adata,color=["tissueProv"],save=Experiment+'_tissueProv',title=Experiment,show=False)

adata.write(Experiment + post_ext)

print(f'{Experiment} Dimension Reduction DONE')



