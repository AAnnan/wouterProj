import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
#import decoupler as dc
import matplotlib.pyplot as plt
#from matplotlib_venn import venn2,venn3
import anndata as ad
import os
from collections import Counter


#Functions
#Hist of the Distribution
def hist_SoloGenesScore(adata,genes,name_genelist,Experiment,raw=False,log_bool=True):
    lay = 'scran_normalization'
    if raw:
        lay = 'No_norm'
    # Create a figure and three subplots
    fig, axes = plt.subplots(1, len(genes), figsize=(len(genes)*5, 5))
    # Plot histograms for each column
    for i, marker_name in enumerate(genes):
        axes[i].hist(list(adata.layers[lay].T[adata.var.index == marker_name].toarray()[0]), bins=100, color='skyblue',log=log_bool) #, ,range=(-0.2, 0.6)
        # Add labels and title
        axes[i].set_xlabel(marker_name + ' Score')
        if log_bool:
            axes[i].set_ylabel('Log Number of cells')
        else:
            axes[i].set_ylabel('Number of cells')
        axes[i].set_title(marker_name + ' Score Distribution')
    plt.tight_layout()
    plt.savefig(f'./figures/{Experiment}_{name_genelist}Score_Hist.pdf')
    plt.close()

#Hist of the Distribution
def hist_gene_sigs(adata,marker_dict,name_dict,thresh_gs,Experiment,log_bool=True):
    # Create a figure and three subplots
    if len(marker_dict) > 3:
        fig, axes = plt.subplots(2, len(marker_dict)//2, figsize=((len(marker_dict)//2)*5, 10))
        # Plot histograms for each column
        pos = 0
        for i in range(axes.shape[0]):
            for j in range(axes.shape[1]):
                
                axes[i, j].hist(adata.obs[list(marker_dict)[pos]], bins=100, color='skyblue',log=log_bool) #, ,range=(-0.2, 0.6)
                if thresh_gs:
                    axes[i, j].axvline(x=thresh_gs, color='red', linestyle='--', linewidth=1)

                # Add labels and title
                axes[i, j].set_xlabel(list(marker_dict)[pos] + ' Score')
                if log_bool:
                    axes[i, j].set_ylabel('Log Number of cells')
                else:
                    axes[i, j].set_ylabel('Number of cells')
                axes[i, j].set_title(list(marker_dict)[pos] + ' Score Distribution')
                pos+=1
    else:
        fig, axes = plt.subplots(1, len(marker_dict), figsize=(len(marker_dict)*5, 5))
            # Plot histograms for each column
        for i, marker_name in enumerate(list(marker_dict)):
            axes[i].hist(adata.obs[marker_name], bins=100, color='skyblue',log=log_bool) #, ,range=(-0.2, 0.6)
            if thresh_gs:
                axes[i].axvline(x=thresh_gs, color='red', linestyle='--', linewidth=1)


            # Add labels and title
            axes[i].set_xlabel(marker_name + ' Score')
            if log_bool:
                axes[i].set_ylabel('Log Number of cells')
            else:
                axes[i].set_ylabel('Number of cells')
            axes[i].set_title(marker_name + ' Score Distribution')

    plt.tight_layout()
    plt.savefig(f'./figures/{Experiment}_{name_dict}Score_Hist.pdf')
    plt.close()
    return 0

def boxplot_gene_sigs(mod,gsetName,modality,Experiment,clusters='leiden_res0.25'):
    # Calculate mean values for each cluster
    mean_values = mod.obs.groupby(clusters,observed=False)[gsetName].median()
    sorted_clusters = mean_values.sort_values(ascending=False).index.tolist()

    # Create a boxplot using seaborn
    sns.set(style="whitegrid")
    sns.boxplot(x=mod.obs[clusters], y=mod.obs[gsetName], order=sorted_clusters, hue=mod.obs[clusters], dodge=False)

    # Add titles and labels
    plt.title(f'Boxplot of {gsetName} Scores by {clusters}\n{Experiment}')
    plt.xlabel(f'Clusters {clusters}')
    plt.ylabel(f'{gsetName} Score')
    plt.xticks(rotation=90)

    plt.tight_layout()
    plt.savefig(f'./figures/{Experiment}_{modality}boxplots_{clusters}_{gsetName}Score.pdf')
    plt.close()

def get_highscore_clusters(adata_obj,gene_obs,leires='leiden_res2.5',medianthresh=0):
    HighScoreClusters = {}
    for clus in adata_obj.obs[leires].cat.categories:
        MedClus = np.median(adata_obj.obs[gene_obs][adata_obj.obs[leires]==clus])
        print(f'Median Score {gene_obs} {clus}: {MedClus}')
        if MedClus > medianthresh:
            HighScoreClusters[clus] = MedClus
    return HighScoreClusters

def Subtypes_Thresh_ScatterL1L2B(adat,name,Experiment,thresholds = [i / 100 for i in range(0, 151)],saveplot=True):
# distribution of uniquely called subtypes depending on genescoring thresholds
    n_cells = adat.n_obs/100
    l1 = []
    l2 = []
    b = []
    ambig = []
    NoTypes = []

    for th in thresholds:
        setL1 = set(adat.obs.index[adat.obs['Epi_Luminal_1'] > th])
        setL2 = set(adat.obs.index[adat.obs['Epi_Luminal_2Psca'] > th])
        setb = set(adat.obs.index[adat.obs['Epi_Basal_1'] > th])

        l1.append(len(setL1 - (setL2 | setb)))
        l2.append(len(setL2 - (setL1 | setb)))
        b.append(len(setb - (setL1 | setL2)))

        ambig.append(len(set(setb & setL1 | setb & setL2 | setL1 & setL2)))
        NoTypes.append(len(set(adat.obs.index) - set(setL1 | setL2 | setb)))

    if saveplot:
        plt.scatter(x=thresholds,y=np.array(l1)/n_cells, alpha=0.8, color='#a6cee3', edgecolor='black', linewidth=0.2, label="L1 only",s=plt.rcParams['lines.markersize'] ** 1.5)
        plt.scatter(x=thresholds,y=np.array(l2)/n_cells, alpha=0.8, color='#b2df8a', edgecolor='black', linewidth=0.2, label="L2 only",s=plt.rcParams['lines.markersize'] ** 1.5)
        plt.scatter(x=thresholds,y=np.array(b)/n_cells, alpha=0.8, color='#8a6be2', edgecolor='black', linewidth=0.2, label="Basal only",s=plt.rcParams['lines.markersize'] ** 1.5)
        plt.scatter(x=thresholds,y=np.array(ambig)/n_cells, alpha=0.8, color='#bdbdbd', edgecolor='black', linewidth=0.2, label="Ambig",s=plt.rcParams['lines.markersize'] ** 1.5)
        plt.scatter(x=thresholds,y=np.array(NoTypes)/n_cells, alpha=0.8, color='#fb9a99', edgecolor='black', linewidth=0.2, label="Untyped",s=plt.rcParams['lines.markersize'] ** 1.5)
        #plt.yscale('log')
        plt.ylabel(f'% of Epithelial cells (w/o SV)')
        plt.xlabel('Score threshold for subtype assignment')
        plt.title('Uniquely assigned Epithelial subtypes')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/{Experiment}_{name}_NCells_UniqueSubtypes_Scatter.pdf')
        plt.close()

    return thresholds[np.argmin(np.array(ambig)+np.array(NoTypes))]

def Subtypes_Thresh_ScatterSV(adat,name,Experiment,thresholds = [i / 100 for i in range(0, 151)]):
# distribution of uniquely called subtypes depending on genescoring thresholds
    n_cells = adat.n_obs/100
    bSV = []
    lSV = []
    ambig = []
    NoTypes = []

    for th in thresholds:
        setbSV = set(adat.obs.index[adat.obs['Epi_Basal_SV'] > th])
        setlSV = set(adat.obs.index[adat.obs['Epi_Luminal_SV'] > th])

        bSV.append(len(setbSV - setlSV))
        lSV.append(len(setlSV - setbSV))

        ambig.append(len(set(setbSV & setlSV)))
        NoTypes.append(len(set(adat.obs.index) - set(setbSV | setlSV)))

    plt.scatter(x=thresholds,y=np.array(bSV)/n_cells, alpha=0.8, color='#a6cee3', edgecolor='black', linewidth=0.2, label="Basal SV only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(lSV)/n_cells, alpha=0.8, color='#b2df8a', edgecolor='black', linewidth=0.2, label="Luminal SV only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(ambig)/n_cells, alpha=0.8, color='#bdbdbd', edgecolor='black', linewidth=0.2, label="Ambig",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(NoTypes)/n_cells, alpha=0.8, color='#fb9a99', edgecolor='black', linewidth=0.2, label="Untyped",s=plt.rcParams['lines.markersize'] ** 1.5)
    #plt.yscale('log')
    plt.ylabel(f'% of Epithelial cells')
    plt.xlabel('Score threshold for subtype assignment')
    plt.title('Uniquely assigned Epithelial subtypes')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figures/{Experiment}_{name}_NCells_SV_Scatter.pdf')
    plt.close()

    return thresholds[np.argmin(np.array(ambig)+np.array(NoTypes))]

#VENN 3 SUBTYPES
def venn3_subtypes(adat,marker_names_list,threshold,Experiment,nameScoring):
    for gset_name in marker_names_list:
        adat.obs[f'{gset_name}_Type'] = adat.obs[gset_name] > threshold

    class_sets = []
    for gset_name in marker_names_list:
        class_sets.append(set(adat.obs.index[adat.obs[f'{gset_name}_Type']]))

    venn_diagram = venn3(class_sets, set_labels=marker_names_list)
    plt.title("Overlap in Epithelial Subtypes")
    plt.tight_layout()
    plt.savefig(f'./figures/{Experiment}_Epi_Subtypes_{nameScoring}_Venn3.pdf')
    plt.close()

def quant_singleSeries(pdCat):
    counter = Counter(pdCat)
    newSeries = pd.Series(pdCat.tolist())
    
    for typ in counter.keys():
        newSeries[newSeries == typ] = f'{typ} {counter[typ]}cells'
    
    return pd.Categorical(newSeries,categories=[f'{typ} {counter[typ]}cells' for typ in counter.keys()])

# L1/L2 Signature
def karthaus_plotting(adat,timepoint,Experimenti,TypeMethod,nameScoring):
    print(f"Plotting L1/L2 {Experimenti} {TypeMethod} {timepoint} {nameScoring}")
    
    minplot = min(np.min(adat.obs['Epi_Luminal_1']),np.min(adat.obs['Epi_Luminal_2Psca']))
    maxplot = max(np.max(adat.obs['Epi_Luminal_1']),np.max(adat.obs['Epi_Luminal_2Psca']))
    plt.xlim(xmin=minplot, xmax=maxplot)
    plt.ylim(ymin=minplot, ymax=maxplot)

    adat = adat[adat.obs['timePoint']==timepoint].copy()

    vmaxl1 = np.percentile(adat.obs['Epi_Luminal_1'], 90)
    vmaxl2 = np.percentile(adat.obs['Epi_Luminal_2Psca'], 90)

    l1x = adat[adat.obs[TypeMethod].str.startswith('L1')].obs['Epi_Luminal_1']
    l1y = adat[adat.obs[TypeMethod].str.startswith('L1')].obs['Epi_Luminal_2Psca']
    l2x = adat[adat.obs[TypeMethod].str.startswith('L2')].obs['Epi_Luminal_1']
    l2y = adat[adat.obs[TypeMethod].str.startswith('L2')].obs['Epi_Luminal_2Psca']

    if 'DiffAmbig' in TypeMethod:
        ambigx = adat[adat.obs[TypeMethod].str.startswith('Ambig_L1|L2')].obs['Epi_Luminal_1']
        ambigy = adat[adat.obs[TypeMethod].str.startswith('Ambig_L1|L2')].obs['Epi_Luminal_2Psca']
    else:
        ambigx = adat[adat.obs[TypeMethod].str.startswith('Ambig')].obs['Epi_Luminal_1']
        ambigy = adat[adat.obs[TypeMethod].str.startswith('Ambig')].obs['Epi_Luminal_2Psca']

    # Create a scatter plot with different colors for x and y values
    plt.scatter(l1x, l1y,c=l1x,cmap='Reds', s=8,label=f'L1 {len(l1x)}cells',vmin=0,vmax=vmaxl1)
    plt.scatter(l2x, l2y,c=l2y,cmap='Blues', s=8,label=f'L2 {len(l2y)}cells',vmin=0,vmax=vmaxl2)
    if 'Diff' in TypeMethod:
        plt.scatter(ambigx, ambigy,c=ambigx+ambigy,cmap='Greys', s=8,alpha=0.8, label=f'Ambig L1|L2 {len(ambigx)}cells')
    else:
        plt.scatter(ambigx, ambigy,c=ambigx+ambigy,cmap='Greys', s=8,alpha=0.8, label=f'Ambig {len(ambigx)}cells')
    
    # Add labels and title
    plt.xlabel('L1 Signature Score')
    plt.ylabel('L2 Signature Score')
    plt.title(f'{Experimenti} {TypeMethod} {timepoint} {nameScoring}')
    plt.legend()
    
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"./figures/L1L2Sig_{Experimenti}_{TypeMethod}_{timepoint}_{nameScoring}.pdf")
    plt.close()
    

def hist_gene_sigs(adata,marker_dict,name_dict,thresh_gs,Experiment,log_bool=True):
        fig, axes = plt.subplots(1, len(marker_dict), figsize=(len(marker_dict)*5, 5))
            # Plot histograms for each column
        for i, marker_name in enumerate(list(marker_dict)):
            axes[i].hist(adata.obs[marker_name], bins=100, color='skyblue',log=log_bool) #, ,range=(-0.2, 0.6)
            if thresh_gs:
                axes[i].axvline(x=thresh_gs, color='red', linestyle='--', linewidth=1)


            # Add labels and title
            axes[i].set_xlabel(marker_name + ' Score')
            if log_bool:
                axes[i].set_ylabel('Log Number of cells')
            else:
                axes[i].set_ylabel('Number of cells')
            axes[i].set_title(marker_name + ' Score Distribution')

        plt.tight_layout()
        plt.savefig(f'./figures/{Experiment}_{name_dict}Score_Hist.pdf')
        plt.close()
        return 0

# def score_genes_with_weights(
# adata: AnnData,
# gene_list: Sequence[str] | pd.Index[str],
# *,
# ctrl_size: int = 50,
# gene_pool: Sequence[str] | pd.Index[str] | None = None,
# n_bins: int = 25,
# score_name: str = "score",
# random_state: AnyRandom = 0,
# copy: bool = False,
# use_raw: bool | None = None,
# ) -> AnnData | None:
# """\
# Score a set of genes [Satija15]_.

# The score is the average expression of a set of genes subtracted with the
# average expression of a reference set of genes. The reference set is
# randomly sampled from the `gene_pool` for each binned expression value.

# This reproduces the approach in Seurat [Satija15]_ and has been implemented
# for Scanpy by Davide Cittaro.

# Parameters
# ----------
# adata
#     The annotated data matrix.
# gene_list
#     The list of gene names used for score calculation.
# ctrl_size
#     Number of reference genes to be sampled from each bin. If `len(gene_list)` is not too
#     low, you can set `ctrl_size=len(gene_list)`.
# gene_pool
#     Genes for sampling the reference set. Default is all genes.
# n_bins
#     Number of expression level bins for sampling.
# score_name
#     Name of the field to be added in `.obs`.
# random_state
#     The random seed for sampling.
# copy
#     Copy `adata` or modify it inplace.
# use_raw
#     Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.

#     .. versionchanged:: 1.4.5
#        Default value changed from `False` to `None`.

# Returns
# -------
# Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following field:

# `adata.obs[score_name]` : :class:`numpy.ndarray` (dtype `float`)
#     Scores of each cell.

# Examples
# --------
# See this `notebook <https://github.com/scverse/scanpy_usage/tree/master/180209_cell_cycle>`__.
# """
# start = logg.info(f"computing score {score_name!r}")
# adata = adata.copy() if copy else adata
# use_raw = _check_use_raw(adata, use_raw)


# if random_state is not None:
#     np.random.seed(random_state)


# var_names = adata.raw.var_names if use_raw else adata.var_names
# gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
# genes_to_ignore = gene_list.difference(var_names, sort=False)  # first get missing
# gene_list = gene_list.intersection(var_names)  # then restrict to present
# if len(genes_to_ignore) > 0:
#     logg.warning(f"genes are not in var_names and ignored: {genes_to_ignore}")
# if len(gene_list) == 0:
#     raise ValueError("No valid genes were passed for scoring.")


# if gene_pool is None:
#     gene_pool = pd.Index(var_names, dtype="string")
# else:
#     gene_pool = pd.Index(gene_pool, dtype="string").intersection(var_names)
# if len(gene_pool) == 0:
#     raise ValueError("No valid genes were passed for reference set.")


# # Trying here to match the Seurat approach in scoring cells.
# # Basically we need to compare genes against random genes in a matched
# # interval of expression.


# _adata = adata.raw if use_raw else adata
# _adata_subset = (
#     _adata[:, gene_pool] if len(gene_pool) < len(_adata.var_names) else _adata
# )
# # average expression of genes
# if issparse(_adata_subset.X):
#     obs_avg = pd.Series(
#         np.array(_sparse_nanmean(_adata_subset.X, axis=0)).flatten(),
#         index=gene_pool,
#     )
# else:
#     obs_avg = pd.Series(np.nanmean(_adata_subset.X, axis=0), index=gene_pool)


# # Sometimes (and I don't know how) missing data may be there, with nansfor
# obs_avg = obs_avg[np.isfinite(obs_avg)]


# n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
# obs_cut = obs_avg.rank(method="min") // n_items
# control_genes = pd.Index([], dtype="string")


# # now pick `ctrl_size` genes from every cut
# for cut in np.unique(obs_cut.loc[gene_list]):
#     r_genes: pd.Index[str] = obs_cut[obs_cut == cut].index
#     if ctrl_size < len(r_genes):
#         r_genes = r_genes.to_series().sample(ctrl_size).index
#     control_genes = control_genes.union(r_genes.difference(gene_list))


# X_list = _adata[:, gene_list].X
# if issparse(X_list):
#     X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()
# else:
#     X_list = np.nanmean(X_list, axis=1, dtype="float64")


# X_control = _adata[:, control_genes].X
# if issparse(X_control):
#     X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
# else:
#     X_control = np.nanmean(X_control, axis=1, dtype="float64")


# score = X_list - X_control


# adata.obs[score_name] = pd.Series(
#     np.array(score).ravel(), index=adata.obs_names, dtype="float64"
# )


# logg.info(
#     "    finished",
#     time=start,
#     deep=(
#         "added\n"
#         f"    {score_name!r}, score of gene set (adata.obs).\n"
#         f"    {len(control_genes)} total control genes are used."
#     ),
# )
# return adata if copy else None
