import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix, issparse
import anndata as ad
import os

#For Correction of ambient RNA
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

DirRNA = '/mnt/etemp/ahrmad/wouter/batch_RNA'
Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
### TO CHANGE IF CHANGING SAMPLES ###
Experiment='Wouter21_ENZ_AP'
### TO CHANGE IF CHANGING SAMPLES ###

qc_ext = '_qc.h5ad'
resDir = '/mnt/etemp/ahrmad/wouter/refs'
Feat_Dim_Red_ext = '_Norm_FeatSel_DimRed.h5ad'

Samples = [d for d in os.listdir(Dir10x) if d.startswith('WK')]
sample_dict = {'Enzymatic Digestion':[],'Singulator':[]}
for sample in Samples:
	if '1350' in sample:
		sample_dict['Enzymatic Digestion'].append(sample)
	else:
		sample_dict['Singulator'].append(sample)
if 'ENZ' in Experiment:
	Samples = sample_dict['Enzymatic Digestion']
elif 'SING' in Experiment:
	Samples = sample_dict['Singulator']
if 'ENZ_AP' in Experiment:
	Samples = [s for s in Samples if 'AP' in s]


print(f'Loading Metadata...')
BC_dict = {}
for sample in Samples:
	atacdf = pd.read_csv(f'{resDir}/csv/{sample}_doublet_scores_ATACPeaks_Self.csv')
	gexdf = pd.read_csv(f'{resDir}/csv/{sample}_doublet_scores_GEX.csv')
	merged_df_all = atacdf.merge(gexdf, on='obs_names', how='inner')
	merged_df_all['doublet_class'] = 'WHATAMI'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Doublet ATAC Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet GEX Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
	merged_df_all_singlet = merged_df_all[merged_df_all['doublet_class'].str.contains('Only')]
	singlets = list(merged_df_all_singlet['obs_names'])
	BC_dict[sample] = singlets

seqDate_dict = {}
with open(f'{resDir}/sequencing_batches.txt', 'r') as seqbatchfile:
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

#Concatenation
adata_list = []
print(f'Loading AnnData...')
for sample in Samples:
    print(f'Reading {sample}')
    a = sc.read_h5ad(os.path.join(DirRNA,sample+qc_ext))
    a.obs['tissueProv'] = tissueProv_dict[sample]
    a.obs['seqDate'] = seqDate_dict[sample]
    a.obs['timePoint'] = timePoint_dict[sample]
    a.obs['isoMeth'] = isoMeth_dict[sample]
    a.obs['mouseID'] = mouseID_dict[sample]
    
    #Filter doublets per sample
    a = a[np.isin(a.obs.index, BC_dict[sample])].copy()
    adata_list.append(a)
    #del a

print(f'Loading Done\nConcatenating...')
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
plt.savefig(Experiment+'_Norm.pdf')
plt.close()

print('DONE with SCRAN Normalization')

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
plt.savefig(Experiment+'_Disp.pdf')
plt.close()

print(f'{Experiment} Feature Selection with Deviance DONE')


adata.X = adata.layers["scran_normalization"]

#PCA
# setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_variance_ratio(adata, log=True,save=Experiment+'_PCA_Variance',show=False)
print(f'{Experiment} Dimension Reduction DONE')

sc.pp.neighbors(adata, use_rep="X_pca",n_neighbors=30, n_pcs=30)
sc.tl.umap(adata)

for resLeiden in [.25,.5,1,1.25,1.5,2,2.5,3]:
	print(f'Leiden clustering at {resLeiden} resolution')
	sc.tl.leiden(adata, key_added=f"leiden_res{resLeiden}", resolution=resLeiden)
	sc.pl.umap(adata,color=f"leiden_res{resLeiden}",legend_loc="on data",save=f"{Experiment}_leiden_res{resLeiden}",title=Experiment,show=False)

#OTHER UMAP PLOTS
sc.pl.umap(adata,color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],save=Experiment+'_UMAP_QC',show=False)
sc.pl.umap(adata,color=["batch"],save=Experiment+'_batch',title=Experiment,show=False)
sc.pl.umap(adata,color=["timePoint"],save=Experiment+'_timePoint',title=Experiment,show=False)
sc.pl.umap(adata,color=["isoMeth"],save=Experiment+'_isoMeth',title=Experiment,show=False)
sc.pl.umap(adata,color=["mouseID"],save=Experiment+'_mouseID',title=Experiment,show=False)
sc.pl.umap(adata,color=["seqDate"],save=Experiment+'_seqDate',title=Experiment,show=False)
sc.pl.umap(adata,color=["tissueProv"],save=Experiment+'_tissueProv',title=Experiment,show=False)

isoMeth2 = []
for el in adata.obs.index:
    if '1350' in el:
        isoMeth2.append('Enzymatic Digestion')
    elif 'Citrate' in el or 'Castrate' in el or 'Contrl' in el or 'Day3' in el:
        isoMeth2.append('Singulator (Putative)')
    elif 'BL6_I' in el or 'BL6_AP' in el:
        isoMeth2.append('Singulator')
adata.obs["isoMeth2"] = isoMeth2
sc.pl.umap(adata,color=["isoMeth2"],save=Experiment+'_isoMethPut',title=Experiment,show=False)
print(f'{Experiment} UMAP DONE')

print('Writing final h5ad')
adata.write(Experiment + Feat_Dim_Red_ext)



print('Annotation')
adata = sc.read_h5ad(Experiment + Feat_Dim_Red_ext)
### Annotation
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import numba
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
sc.set_figure_params(figsize=(5, 5))

Feat_Dim_Red_ext = '_Norm_FeatSel_DimRed.h5ad'
inputNormH5ad = [d.removesuffix(Feat_Dim_Red_ext) for d in os.listdir('.') if d.endswith(Feat_Dim_Red_ext)]

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

#UMAP GENE COLORS
adata.var.index = adata.var.index.str.upper()

markerGenes_in_data = dict()
for ct, markers in markerGenes.items():
    markers_found = list()
    for marker in markers:
        if marker.upper() in adata.var.index:
            markers_found.append(marker.upper())
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
    plt.savefig(Experiment+f'_UMAP_{ct.upper()}.pdf')
    plt.close()

#DOTPLOT
sc.pl.dotplot(
    adata,
    groupby="leiden_res0.25",
    var_names=markerGenes_in_data,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)
plt.suptitle(Experiment)
plt.subplots_adjust(bottom = 0.25)
plt.tight_layout()
plt.savefig(Experiment+f'_DOTPLOT_LUMINAL.pdf')
plt.close()

#Celltyping based on Paper genes 
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L1']], ctrl_size=50, n_bins=25, score_name='L1_genes_Paper')
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L2']], ctrl_size=50, n_bins=25, score_name='L2_genes_Paper')
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L3']], ctrl_size=50, n_bins=25, score_name='L3_genes_Paper')
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L']], ctrl_size=50, n_bins=25, score_name='L_genes_Paper')

sc.pl.umap(adata,color=["L1_genes_Paper"],save=Experiment+'_L1_genes_fromMain',title=f'{Experiment}_L1 sig paper',show=False)
sc.pl.umap(adata,color=["L2_genes_Paper"],save=Experiment+'_L2_genes_fromMain',title=f'{Experiment}_L2 sig paper',show=False)
sc.pl.umap(adata,color=["L3_genes_Paper"],save=Experiment+'_L3_genes_fromMain',title=f'{Experiment}_L3 sig paper',show=False)
sc.pl.umap(adata,color=["L_genes_Paper"],save=Experiment+'_L_genes_fromMain',title=f'{Experiment}_L sig paper',show=False)

#Celltyping based on Supp gene lists
gene_set_table = f'{resDir}/table_s8_summary.txt'
gene_set_df = pd.read_csv(gene_set_table, sep='\t')

for gene_set in gene_set_df.columns:
	ctrl = 50
	print(f'{gene_set} S8')
	genes = gene_set_df[gene_set].str.upper()
	sc.tl.score_genes(adata, genes, ctrl_size=ctrl, n_bins=25, score_name=gene_set)
	sc.pl.umap(adata,color=gene_set,save=f'{Experiment}_{gene_set}_SuppS8_ctrl{ctrl}',title=f'{Experiment}_{gene_set}',show=False)



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
	genes = pd.Series(sub_ct_marker_list[gene_set]).str.upper()
	sc.tl.score_genes(adata, genes, ctrl_size=50, n_bins=25, score_name=gene_set)
	sc.pl.umap(adata,color=gene_set,save=f'{Experiment}_{gene_set}_sub_ct_marker',title=f'{Experiment}_{gene_set}_sub_ct_marker',show=False)

for gene_set in main_ct_marker_list.keys():
	print(f'{gene_set} DanDict')
	genes = pd.Series(main_ct_marker_list[gene_set]).str.upper()
	sc.tl.score_genes(adata, genes, ctrl_size=50, n_bins=25, score_name=gene_set)
	sc.pl.umap(adata,color=gene_set,save=f'{Experiment}_{gene_set}_main_ct_marker',title=f'{Experiment}_{gene_set}_main_ct_marker',show=False)

#L1_genes = typing_gene_df['Epi_Luminal_1'].str.upper()
#L2_genes = typing_gene_df['Epi_Luminal_2Psca'].str.upper()
#L3_genes = typing_gene_df['Epi_Luminal_3Foxi1'].str.upper()
#L_genes = typing_gene_df['Epi_Luminal_SV'].str.upper() #seminal vesicule
#
#sc.tl.score_genes(adata, L1_genes, ctrl_size=200, n_bins=25, score_name='L1_genes')
#sc.tl.score_genes(adata, L2_genes, ctrl_size=200, n_bins=25, score_name='L2_genes')
#sc.tl.score_genes(adata, L3_genes, ctrl_size=200, n_bins=25, score_name='L3_genes')
#sc.tl.score_genes(adata, L_genes, ctrl_size=200, n_bins=25, score_name='L_genes')
#
# of genes in lists not found in adata
#sum(~np.isin(pd.concat([L_genes,L1_genes,L2_genes,L3_genes]),adata.var.index))
#
#sc.pl.umap(adata,color=["L1_genes"],save=Experiment+'_L1_genes_fromSupp',title=f'{Experiment}_L1 sig',show=False)
#sc.pl.umap(adata,color=["L2_genes"],save=Experiment+'_L2_genes_fromSupp',title=f'{Experiment}_L2 sig',show=False)
#sc.pl.umap(adata,color=["L3_genes"],save=Experiment+'_L3_genes_fromSupp',title=f'{Experiment}_L3 sig',show=False)
#sc.pl.umap(adata,color=["L_genes"],save=Experiment+'_L_genes_fromSupp',title=f'{Experiment}_L sig',show=False)


print('Writing final h5ad')
adata.write(Experiment + Feat_Dim_Red_ext)
sc.pl.highest_expr_genes(adata, n_top=10,save=Experiment, show=False)
#adata = sc.read_h5ad(Experiment + Feat_Dim_Red_ext)


#import snapatac2 as snap
#snap.pl.umap(adata,color="batch",out_file=Experiment+'_batch.html',interactive=False)

#main_ct_marker_list = list(Epithelial = c("Epcam"), Immune = "Ptprc") # present
# sub_ct_marker_list = list(Bcells = c("Cd19","Ms4a1"), Tcells = c("Cd3g","Cd3d","Cd3e","Cd247"), macrophages = c("Cd14", "Aif1"), endothelium = c("Pecam1","Vwf"),
#       lymphatic_endothelium = c("Pecam1","Prox1" ), glia = c("Sox10"), myofibroblast = c("Acta2","Myh11","Rspo3" ), smooth_muscle = c("Acta2", "Notch3"),
#       mesenchymal_1 = c("Col5a2", "Lama2", "Zeb1", "Wnt2", "Wnt6", "Wnt10a", "Rorb"), mesenschymal_2 = c("Col5a2", "Lama2", "Zeb1", "Sult1e1", "Fgf10", "Rspo1"),
#       basal = c("Epcam", "Krt5", "Krt14", "Trp63"), seminal_vesicle_basal = c("Epcam", "Pax2", "Krt5", "Krt14","Trp63", "Calml3"), 
#       luminal_1 = c("Epcam", "Cd24a", "Krt8", "Krt18", "Nkx3-1", "Pbsn", PROM1, CD26/DPP4)
#, luminal_2 = c("Epcam", "Cd24a", "Krt8", "Krt18", "Psca", "Krt4", "Tacstd2", "Ly6a"),
#       luminal_3 = c("Epcam", "Cd24a", "Krt8", "Krt18", "Foxi1", "Atp6v1g3", "Atp6v1b1"),
# seminal_vesicle_luminal = c( "Epcam", "Pax2", "Krt8", "Krt18", "Pate4"),
#        seminal_vesicle_ionocyte = c("Epcam", "Pax2", "Foxi1") ) # already alias-corrected
### big absences from 'measured_genes_all': none
### big absences from 'measured_genes' (Experiment-specific): Cd19, Cd3*, Cd31/Pecam1, Sox10, some Wnt, Cd24a, ("Foxi1", "Atp6v1g3", "Atp6v1b1") discriminating the luminal_3, 





