import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix, issparse
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

Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
PathMat = '/outs/filtered_feature_bc_matrix.h5'
Samples = [d for d in os.listdir(Dir10x) if d.startswith('WK')]

### Quality Control 
pct_mt_count_thresh = 8
mt_MADS_thresh = 3
count_thresh = 4
filt_top20genes_outliers = True

corr_ambient_rna = 'CellBender'
#SoupX: https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#correction-of-ambient-rna
#CellBender: https://cellbender.readthedocs.io/en/latest/introduction/index.html

qc_ext = '_qc.h5ad'
logFile = 'Wouter_scRNA_log_QC'

for sample in Samples:
    print(f"{sample}")
    adata = sc.read_h5ad(sample + qc_ext)
    sc.pl.violin(adata, 'pct_counts_mt',save=sample+'-MT',title=sample,show=False)
    sc.pl.violin(adata, 'pct_counts_ribo',save=sample+'-Ribo',title=sample,show=False)
    sc.pl.violin(adata, 'pct_counts_hb',save=sample+'-HB',title=sample,show=False)

l = open(logFile, 'a')
for sample in Samples:
    print(f'\n###################\n Processing {sample}...\n###################')
    l.write(f'\n###################\n Processing {sample}...\n###################')
    mat10x = Dir10x + sample + PathMat
    adata = sc.read_10x_h5(filename=mat10x)

    adata.var_names_make_unique()
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.startswith("Hb")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20])
    #We have replaced annotations for lincRNAs Gm42418 and
    #AY036118 with a single contiguous gene annotation for the rRNA element Rn45s. This locus harbors an Rn45s repeat
    #as reflected in RefSeq, such that contaminating 18S rRNA in our library preparations may lead to inflated expression
    #counts for these lincRNAs
    sc.pl.highest_expr_genes(adata, n_top=15,save=sample+'_HighExpGenes', show=False)
    sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    plt.title(sample)
    plt.savefig(sample+'_TotCounts.pdf')
    plt.close()

    #sc.pl.violin(adata, 'pct_counts_mt',save=sample+'-MT',title=sample,show=False)
    #sc.pl.violin(adata, 'pct_counts_ribo',save=sample+'-Ribo',title=sample,show=False)
    #sc.pl.violin(adata, 'pct_counts_hb',save=sample+'-HB',title=sample,show=False)
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",save=sample+'_QC_mt',title=sample,show=False)
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_ribo",save=sample+'_QC_ribo',title=sample,show=False)
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_hb",save=sample+'_QC_hb',title=sample,show=False)

    #### Filtering low quality reads
    #with automatic thresholding based on MAD (median absolute deviations)
    #sc-best-practices.org/preprocessing_visualization/quality_control.html
    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    #log1p_total_counts, log1p_n_genes_by_counts and pct_counts_in_top_20_genes filtered with a threshold of 5 MADs
    adata.obs["outlier_total_c"] = (is_outlier(adata, "log1p_total_counts", count_thresh))
    adata.obs["outlier_n_genesC"] = (is_outlier(adata, "log1p_n_genes_by_counts", count_thresh))
    outlier_total_c = np.sum(adata.obs["outlier_total_c"])
    outlier_n_genesC = np.sum(adata.obs["outlier_n_genesC"])

    if filt_top20genes_outliers:
        adata.obs["outlier_top_20_genes"] = (is_outlier(adata, "pct_counts_in_top_20_genes", count_thresh))
        outlier_top_20_genes = np.sum(adata.obs["outlier_top_20_genes"])
    else:
        adata.obs["outlier_top_20_genes"] = False
        outlier_top_20_genes = np.sum((is_outlier(adata, "pct_counts_in_top_20_genes", count_thresh)))

    #pct_counts_Mt is filtered with mt_MADS_thresh MADs
    #cells with a percentage of mitochondrial counts exceeding 15% are filtered out
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_MADS_thresh) & (adata.obs["pct_counts_mt"] > pct_mt_count_thresh)

    mt_outlier = np.sum(adata.obs["mt_outlier"])

    l.write(f"\nTotal number of cells: {adata.n_obs}")
    l.write(f"\nRemoved {outlier_total_c} cells (log1p_total_counts outliers: {count_thresh} MADS)")
    l.write(f"\nRemoved {outlier_n_genesC} cells (log1p_n_genes_by_counts outliers: {count_thresh} MADS)")
    l.write(f'\nFound {outlier_top_20_genes} pct_counts_in_top_20_genes outlier cells ({count_thresh} MADS). Filtering is {"ON" if filt_top20genes_outliers else "OFF"}')
    l.write(f"\nRemoved {mt_outlier} cells (mt outliers: pct_counts_mt > {mt_MADS_thresh} MADS & pct_mt_count > {pct_mt_count_thresh}%)") 
    cellLowQ = adata.n_obs
    
    #ACTUAL FILTERING
    #filter AnnData object
    adata = adata[(~adata.obs.outlier_total_c) & (~adata.obs.outlier_n_genesC) & (~adata.obs.outlier_top_20_genes) & (~adata.obs.mt_outlier)].copy()
    
    l.write(f"\nNumber of cells after filtering of low quality cells: {adata.n_obs}. ")
    cellHighQ = adata.n_obs
    Cellrm = f'{((cellLowQ-cellHighQ)/cellLowQ) * 100:.2f}% ({(cellLowQ-cellHighQ)})'
    l.write(f'\n{Cellrm} low quality cells filtered out')

    ####Correction of ambient RNA
#    if corr_ambient_rna == 'CellBender':
#        cellbender remove-background \
#             --input mat10x.replace("filtered_", "raw_") \
#             --output f"{sample}_CellBender.h5" \
#             --expected-cells adata.n_obs \
#             --total-droplets-included 6000

    if corr_ambient_rna == 'SoupX':
        #SoupX works better on clustered data (even if loose clusters)
        adata_pp = adata.copy()
        sc.pp.normalize_per_cell(adata_pp)
        sc.pp.log1p(adata_pp)
        sc.pp.pca(adata_pp)
        sc.pp.neighbors(adata_pp)
        sc.tl.leiden(adata_pp, key_added="soupx_groups")

        # Preprocess variables for SoupX
        soupx_groups = adata_pp.obs["soupx_groups"]
        del adata_pp

        cells = adata.obs_names
        genes = adata.var_names
        data = adata.X.T
        adata_raw = sc.read_10x_h5(filename=mat10x.replace("filtered_", "raw_"))
        adata_raw.var_names_make_unique()
        data_tod = adata_raw.X.T
        del adata_raw

        %R -i data -i data_tod -i genes -i cells -i soupx_groups
        %R library(SoupX)
        # specify row and column names of data
        %R rownames(data) = genes
        %R colnames(data) = cells
        # ensure correct sparse format for table of counts and table of droplets
        %R data <- as(data, "sparseMatrix")
        %R data_tod <- as(data_tod, "sparseMatrix")
        # Generate SoupChannel Object for SoupX 
        %R sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)
        # Add extra meta data to the SoupChannel object
        %R soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
        %R sc = setSoupProfile(sc, soupProf)
        # Set cluster information in SoupChannel
        %R sc = setClusters(sc, soupx_groups)
        # Estimate contamination fraction
        %R sc  = autoEstCont(sc, doPlot=FALSE)
        # Infer corrected table of counts and rount to integer
        %R outR = adjustCounts(sc, roundToInt = TRUE)

        adata.layers["counts"] = adata.X
        out = %Rget outR
        adata.layers["soupX_counts"] = out.T
        adata.X = adata.layers["soupX_counts"]

    ## LET'S NOT FILTER OUT GENES (PROBLEM WITH MERGING LATER)
    # filter out genes that are not detected in at least 10 cells as these are not informative
    # Min 20 cells - filters out 0 count genes
    #l.write(f"\nTotal number of genes: {adata.n_vars}")
    #geneLowQ = adata.n_vars
    #sc.pp.filter_genes(adata, min_cells=10)
    #l.write(f"\nNumber of genes after cell filter: {adata.n_vars}. ")
    #geneHighQ = adata.n_vars
    #Generm = geneLowQ-geneHighQ
    #l.write(f'\n{Generm}/{geneLowQ} genes filtered out (not detected in at least 10 cells)')

    ##Doublet 
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
    singdoub = adata.obs.scDblFinder_class.value_counts()
    Doubrm = f"{(singdoub['doublet'] / (singdoub['doublet'] + singdoub['singlet'])) * 100:.2f}% ({singdoub['doublet']})"
    l.write(f'\n{Doubrm} doublets')

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="scDblFinder_score",save=sample+'_QCFILT_dbScore',title=sample,show=False)
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="scDblFinder_class",save=sample+'_QCFILT_dbClass',title=sample,show=False)

    #Doublet CSV
    df = pd.DataFrame({'sample': sample, 'obs_names': adata.obs_names, 'doublet_score': adata.obs['scDblFinder_score'], 'doublet_class': adata.obs['scDblFinder_class']})
    df.to_csv(f'{sample}_doublet_scores_GEX.csv', index=False)

    adata.write(sample+qc_ext)

    print(f"QC DONE for {sample}")



### SCRAN Normalization
#log1p with Scran estimated size factors
### Feature selection and Dimensionality Reduction
inputH5ad = [d.removesuffix(qc_ext) for d in os.listdir('.') if d.endswith(qc_ext)]

print(f'Loading Doublet Information...')
BC_dict = {}
for sample in Samples:
    atacdf = pd.read_csv(f'./csv/{sample}_doublet_scores_ATACPeaks_Self.csv')
    gexdf = pd.read_csv(f'./csv/{sample}_doublet_scores_GEX.csv')
    merged_df_all = atacdf.merge(gexdf, on='obs_names', how='inner')
    merged_df_all['doublet_class'] = 'WHATAMI'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Doublet ATAC Only'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet GEX Only'
    merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
    
    # Retain QC passing cells that were called singlets by at least 1 modality
    merged_df_all_singlet = merged_df_all[merged_df_all['doublet_class'].str.contains('Only')]
    singlets = list(merged_df_all_singlet['obs_names'])
    BC_dict[sample] = singlets

%R library(scran)
%R library(BiocParallel)
%R library(scry)
xTopDeviantGenes = 5000
Feat_Dim_Red_ext = '_Norm_FeatSel_DimRed.h5ad'

for sample in inputH5ad:
    print(f"Norm + FeatSel + DimRed - {sample}")
    adata = sc.read_h5ad(sample + qc_ext)

    #FILTER doublets
    #adata = adata[adata.obs['scDblFinder_class']=='singlet'].copy()
    adata = adata[adata.obs.index.isin(BC_dict[sample])].copy()

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
    fig.suptitle(sample)
    plt.savefig(sample+'_Norm.pdf')
    plt.close()

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
    plt.title(sample)
    plt.savefig(sample+'_Disp.pdf')
    plt.close()

    print(f'{sample} Feature Selection DONE')

    adata.X = adata.layers["scran_normalization"]

    #PCA
    # setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    sc.pl.pca_variance_ratio(adata, log=True,save=sample+'_PCA_Variance',show=False)
    
    sc.pp.neighbors(adata, use_rep="X_pca",n_neighbors=30, n_pcs=30)
    sc.tl.umap(adata)
    
    sc.pl.umap(adata,color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],save=sample+'_UMAP_QC',title=sample,show=False)

    sc.tl.leiden(adata)
    sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

    sc.pl.umap(adata,color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],legend_loc="on data",save=sample+'_Leiden_Res',title=sample,show=False)
    sc.pl.umap(adata,color=["leiden_res0_25"],legend_loc="on data",save=sample+'_Leiden025',title=sample,show=False)

    adata.write(sample + Feat_Dim_Red_ext)

    print(f'{sample} Dimension Reduction DONE')

### Annotation
import celltypist
from celltypist import models
import scarches as sca
import urllib.request

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
    "L1": ["Pbsn","Nkx3.1","CD26","Dpp4","CD59a","CD133","Prom1"],  #L1
    "L2": ["Sca1","Ly6a","Tacstd2","Trop2","Psca","Krt4","Claudin10"], #L2
    "L3": ["Foxi1","Atp6v1g3","Atp6b1b"], #L3
    "L": ["CD24a","Krt8","Krt18"] #All Luminal cells
    }
#SUPP p.12
#luminal 1, Epcam, CD24a, Krt8, Krt18, Nkx3.1, Pbsn high
#luminal 2, Epcam, CD24a, Krt8, Krt18, Psca, Krt4, Tacst2, Ly6a
#luminal 3 (ionocyte), Epcam, CD24a, Krt8, Krt18, Foxi1, Atp6v1g3, Atp6b1b

sample="WK-1350_R3-2_DP"

#UMAP GENE COLORS
for sample in inputNormH5ad:
    adata = sc.read_h5ad(sample + Feat_Dim_Red_ext)
    adata.var.index = adata.var.index.str.upper()

    markerGenes_in_data = dict()
    for ct, markers in markerGenes.items():
        markers_found = list()
        for marker in markers:
            if marker.upper() in adata.var.index:
                markers_found.append(marker.upper())
        markerGenes_in_data[ct] = markers_found

    for ct in list(markerGenes_in_data.keys()):
        print(f"{sample} - {ct.upper()}")  # print cell subtype name
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

#for i in {1..3}; do
#    pdfjam --nup 2x6 --landscape --a4paper cellTyping/*UMAP*L${i}.pdf --outfile Wouter_ATAC_GeneSig_L${i}.pdf
#done
#pdfjam --nup 2x6 --landscape --a4paper cellTyping/*UMAP*L.pdf --outfile Wouter_ATAC_GeneSig_L.pdf

#DOTPLOT
for sample in inputNormH5ad:
    print(f"{sample}") 
    adata = sc.read_h5ad(sample + Feat_Dim_Red_ext)
    adata.var.index = adata.var.index.str.upper()

    markerGenes_in_data = dict()
    for ct, markers in markerGenes.items():
        markers_found = list()
        for marker in markers:
            if marker.upper() in adata.var.index:
                markers_found.append(marker.upper())
        markerGenes_in_data[ct] = markers_found

    sc.pl.dotplot(
        adata,
        groupby="leiden_res0_25",
        var_names=markerGenes_in_data,
        standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    )
    plt.suptitle(sample)
    plt.subplots_adjust(bottom = 0.25)
    plt.tight_layout()
    plt.savefig(sample+f'_DOTPLOT_LUMINAL.pdf')
    plt.close()








    

