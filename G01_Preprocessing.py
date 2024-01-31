import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix, issparse
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#For R code
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

CellBender = True

if CellBender:
	Dir10x = '/mnt/etemp/ahrmad/wouter/CellBender/'
	PathMat = '_CellBender_filtered.h5'
	cantCellBend = ['WK-1580_BL6_AP_RegenDay1']
else:
	Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
	PathMat = '/outs/filtered_feature_bc_matrix.h5'

Samples = [d.split(PathMat)[0] for d in os.listdir(Dir10x) if d.startswith('WK') if d.endswith(PathMat)]
Samples = Samples + cantCellBend if CellBender else Samples

refDir = '/mnt/etemp/ahrmad/wouter/refs'
if not os.path.exists(refDir+'/csv'):
    os.makedirs(refDir+'/csv')
    print(f"Directory {refDir+'/csv'} created.")
else:
    print(f"Directory {refDir+'/csv'} already exists.")

### Quality Control 
pct_mt_count_thresh = 10
mt_MADS_thresh = 2
count_MADS_thresh = 5
max_HB_Ribo_pct = 70 #%
max_Top20_pct = 85 #%

qc_ext = '_qc.h5ad'
gex_doub_csv_ext = '_doublet_scores_CB_GEX.csv'
logFile = 'Wouter_scRNA_log_QC'

l = open(logFile, 'a')
for sample in Samples:
	print(f'\n###################\n Processing {sample}...\n###################')
	if sample in cantCellBend:
		mat10x = f'/mnt/ndata/daniele/wouter/Processed/CellRangerArc/{sample}/outs/filtered_feature_bc_matrix.h5'
	else:
		mat10x = Dir10x + sample + PathMat
	adata = sc.read_10x_h5(filename=mat10x)
	adata.X = adata.X.astype(np.float64)
	adata.var_names_make_unique()

	###QC METRICS
	# mitochondrial genes |ribosomal genes | hemoglobin genes.
	adata.var["mt"] = adata.var_names.str.startswith("mt-")
	adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
	adata.var["hb"] = adata.var_names.str.startswith("Hb")
	adata.obs['batch']=sample

	sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
	# Prefiltering of cells > 1000 counts and >100 genes with positive counts
	adata = adata[(adata.obs["total_counts"] > 1000) | (adata.obs["n_genes_by_counts"] > 100)].copy()

	# Barplot of the amounts of reads in 15 top expressed genes
	sc.pl.highest_expr_genes(adata, n_top=15,save=sample)
	# Histogram of otal counts in cells 
	sns.displot(adata.obs["total_counts"], bins=100, kde=False)
	plt.title(sample)
	plt.savefig(os.getcwd()+'/figures/'+sample+'_TotCounts.pdf')
	plt.close()

	# Violin plots
	sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_in_top_20_genes',],multi_panel=True,groupby='batch',rotation=0.0000001,save=sample+'_counts_top20.pdf')
	sc.pl.violin(adata, ['pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'],multi_panel=True,groupby='batch',rotation=0.0000001,save=sample+'_MT_Ribo_HB.pdf')

	# Scatter plots of MT, ribo, hb, top_20_genes counts in each cell
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_mt',title=sample+'\npct_counts_mt',save=sample+'_MT.pdf')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_ribo',title=sample+'\npct_counts_ribo',save=sample+'_Ribo.pdf')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_hb',title=sample+'\npct_counts_hb',save=sample+'_HB.pdf')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_in_top_20_genes',title=sample+'\npct_counts_in_top_20_genes',save=sample+'_top_20_genes.pdf')

	###FILTERING LOW QUALITY CELLS
	#with automatic thresholding based on MAD (median absolute deviations)
	#sc-best-practices.org/preprocessing_visualization/quality_control.html
	def is_outlier(adata, metric: str, nmads: int):
	    M = adata.obs[metric]
	    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
	        np.median(M) + nmads * median_abs_deviation(M) < M
	    )
	    return outlier

	#log1p_total_counts, log1p_n_genes_by_counts filtered with a threshold of "count_MADS_thresh" MADs
	adata.obs["outlier_total_c"] = (is_outlier(adata, "log1p_total_counts", count_MADS_thresh))
	adata.obs["outlier_n_genesC"] = (is_outlier(adata, "log1p_n_genes_by_counts", count_MADS_thresh))
	
	# pct_counts_ribo,pct_counts_hb filtered with a hard threshold of "max_HB_Ribo_pct"%
	adata.obs["outlier_Ribo_HB"] = (adata.obs['pct_counts_ribo'] > max_HB_Ribo_pct) | (adata.obs['pct_counts_hb'] >= max_HB_Ribo_pct)
	# pct_counts_in_top_20_genes filtered with a hard threshold of "max_Top20_pct"%
	adata.obs["outlier_top_20_genes"] = adata.obs['pct_counts_in_top_20_genes'] >= max_Top20_pct
	#adata.obs["outlier_top_20_genes"] = (is_outlier(adata, "pct_counts_in_top_20_genes", count_MADS_thresh))

	#pct_counts_Mt is filtered with "mt_MADS_thresh" MADs AND percentage of mitochondrial counts exceeding "pct_mt_count_thresh"
	adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_MADS_thresh) & (adata.obs["pct_counts_mt"] > pct_mt_count_thresh)
	#cells with a percentage of mitochondrial counts exceeding "pct_mt_count_thresh" are labeled
	adata.obs["mt_above_8pct"] = (adata.obs["pct_counts_mt"] > pct_mt_count_thresh)


	ncell_outlier_top_20_genes = np.sum(adata.obs["outlier_top_20_genes"])
	outlier_total_c = np.sum(adata.obs["outlier_total_c"])
	outlier_n_genesC = np.sum(adata.obs["outlier_n_genesC"])
	ncell_mt_outlier_log = np.sum(adata.obs["mt_outlier"])
	ncell_mt_above_8pct_log = np.sum(adata.obs["mt_above_8pct"])
	ncellLowQ_log = adata.n_obs

	#ACTUAL FILTERING
	#filter AnnData object
	adata = adata[(~adata.obs.outlier_total_c) & (~adata.obs.outlier_n_genesC) & (~adata.obs.mt_outlier) & (~adata.obs.outlier_top_20_genes) & (~adata.obs.outlier_Ribo_HB)].copy()
	
	# Log # cells Postfiltering
	#record n cells after filtering
	ncellHighQ_log = adata.n_obs
	#record the diff
	Cellrm_log = f'{((ncellLowQ_log-ncellHighQ_log)/ncellLowQ_log) * 100:.2f}% ({(ncellLowQ_log-ncellHighQ_log)})'

	##Doublet Analysis
	%R library(Seurat)
	%R library(scater)
	%R library(scDblFinder)
	%R library(BiocParallel)
	data_mat = adata.X.T
	%R -i data_mat
	%R set.seed(1)
	%R sce = scDblFinder(SingleCellExperiment(list(counts=data_mat)))
	%R doublet_scoreR = sce$scDblFinder.score
	%R doublet_classR = sce$scDblFinder.class

	adata.obs["scDblFinder_score"] = %Rget doublet_scoreR
	adata.obs["scDblFinder_class"] = %Rget doublet_classR

	#Doublet CSV
	df = pd.DataFrame({'sample': sample, 'obs_names': adata.obs_names, 'doublet_score': adata.obs['scDblFinder_score'], 'doublet_class': adata.obs['scDblFinder_class']})
	df.to_csv(f'{refDir}/csv/{sample}{gex_doub_csv_ext}', index=False)

	# Scatter plots of MT, ribo, hb, top_20_genes counts in each cell post filtering
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_mt',title=sample+'\npct_counts_mt_PostQC',save=sample+'_MT_PostQC.pdf')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_ribo',title=sample+'\npct_counts_ribo_PostQC',save=sample+'_Ribo_PostQC.pdf')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_hb',title=sample+'\npct_counts_hb_PostQC',save=sample+'_HB_PostQC.pdf')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",color='pct_counts_in_top_20_genes',title=sample+'\npct_counts_in_top_20_genes_PostQC',save=sample+'_top_20_genes_PostQC.pdf')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="scDblFinder_score",title=sample+'\n_dbScore_PostQC',save=sample+'_dbScore_PostQC')
	sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="scDblFinder_class",title=sample+'\n_dbClass_PostQC',save=sample+'_dbClass_PostQC')

	# Log doublet
	singdoub = adata.obs.scDblFinder_class.value_counts()
	Doubrm_log = f"{(singdoub['doublet'] / (singdoub['doublet'] + singdoub['singlet'])) * 100:.2f}% ({singdoub['doublet']})"

	adata.write(sample+qc_ext)
	print(f"QC DONE for {sample}")

    l.write(f'\n###################\n Processing {sample}...\n###################')
    l.write(f"\nTotal number of cells: {adata.n_obs}")
    l.write(f"\nRemoved {outlier_total_c} cells (log1p_total_counts outliers: {count_MADS_thresh} MADS)")
    l.write(f"\nRemoved {outlier_n_genesC} cells (log1p_n_genes_by_counts outliers: {count_MADS_thresh} MADS)")
    l.write(f'\nRemoved {ncell_outlier_top_20_genes} pct_counts_in_top_20_genes outlier cells ({count_MADS_thresh} MADS).')
    l.write(f"\nRemoved {ncell_mt_outlier_log} cells (mt outliers: pct_counts_mt > {mt_MADS_thresh} MADS & pct_mt_count > {pct_mt_count_thresh}%)") 
    l.write(f"\nLabeled {ncell_mt_above_8pct_log} cells as high mt cells (pct_mt_count > {pct_mt_count_thresh}%)") 
    l.write(f"\nNumber of cells after filtering of low quality cells: {ncellHighQ_log}. ")
    l.write(f'\n{Cellrm_log} low quality cells filtered out')
    l.write(f'\n{Doubrm_log} doublets')
l.close()



