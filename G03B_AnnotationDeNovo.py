import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import xlsxwriter
#import decoupler as dc
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
from mpl_toolkits.axes_grid1 import make_axes_locatable
import anndata as ad
import os
import warnings
from collections import Counter
from scipy.stats import median_abs_deviation

from G03_Annotation_Functions import *

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'.*")

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
print("Scanpy version:",sc.__version__)

### TO CHECK WHEN CHANGING SAMPLES ###
Experiment='Wouter21_SING_CB'
Experiment='Wouter21_ENZ_CB'
### TO CHECK WHEN CHANGING SAMPLES ###

refDir = '/mnt/etemp/ahrmad/wouter/refs'
#refDir = './refs'

if not os.path.exists(refDir):
    raise ValueError('Check ref folder paths')
gene_set_table = f'{refDir}/table_s8_summary.txt'
cell_cycle_gene_table = f'{refDir}/Mouse_G2M_S_genes_Tirosh2016.csv'

#Epi subtypes Supp gene lists
gset_df = pd.read_csv(gene_set_table, sep='\t')
cc_df = pd.read_csv(cell_cycle_gene_table)

post_ext = '_post.h5ad'

normToUse = 'SCT_data' #'scran_normalization', 'SCT_data', 'SCT_counts', 'log1p_norm', 'No_norm'

#adata import
adata = sc.read_h5ad(Experiment + post_ext)

u_ext = 'PreCluster'
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'_Counts_QC_'+u_ext,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["batch"],save=Experiment+'_batch_'+u_ext,title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'_timePoint_'+u_ext,title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'_tissueProv_'+u_ext,title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'_Counts_QC_'+u_ext,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],save=Experiment+'_PctCounts_QC_'+u_ext,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["scDblFinder_score","doublet_class"],save=Experiment+'_DbltClass_QC_'+u_ext,show=False)

# Keep only intact
if 'SING' in Experiment:
    adata = adata[np.isin(adata.obs['batch'], np.array(['WK-1501_BL6_INTACT_AP_Test3','WK-1585_INTACT_AP_BL6_Contrl','WK-1585_INTACT_AP_BL6_Citrate']))]
if 'ENZ' in Experiment:
    adata = adata[np.isin(adata.obs['batch'], np.array(['WK-1350_I-2_AP','WK-1350_I-1_AP']))]


# Delete keys from adata.obs containing 'leiden'
leiden_keys_obs = [key for key in adata.obs.keys() if 'leiden' in key]
for key in leiden_keys_obs:
    del adata.obs[key]

# Delete keys from adata.uns containing 'leiden'
leiden_keys_uns = [key for key in adata.uns.keys() if 'leiden' in key]
for key in leiden_keys_uns:
    del adata.uns[key]

sc.pp.highly_variable_genes(adata, layer=normToUse)
adata.X = adata.layers[normToUse]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_variance_ratio(adata, log=True,save=Experiment+'_Epi'+'_PCA_Variance',show=False)
sc.pp.neighbors(adata, use_rep="X_pca",n_neighbors=15, n_pcs=14)
sc.tl.umap(adata)

# redo umaps
u_ext = 'PostCluster'
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'_Counts_QC_'+u_ext,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["batch"],save=Experiment+'_batch_'+u_ext,title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'_timePoint_'+u_ext,title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'_tissueProv_'+u_ext,title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'_Counts_QC_'+u_ext,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],save=Experiment+'_PctCounts_QC_'+u_ext,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["scDblFinder_score","doublet_class"],save=Experiment+'_DbltClass_QC_'+u_ext,show=False)


leiRs = [0.25,1]

for leiR in leiRs:
    leiRes = f'leiden_res{leiR}'
    sc.tl.leiden(adata, key_added=leiRes, resolution=leiR)
    sc.pl.umap(adata,color=leiRes,legend_loc="on data",save=f"{Experiment}_{leiRes}_{u_ext}",title=Experiment,show=False)


for leiR in leiRs:
    sc.tl.rank_genes_groups(adata, groupby=leiRes, method="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(adata, groupby=leiRes, standard_scale="var", n_genes=5, swap_axes=True,show=False,save=Experiment+f'_Intact_DGE_{leiRes}.pdf')
    #sc.get.rank_genes_groups_df(adata, group="2").head(10)

    gset_dict = dict()
    for gset_name in gset_df.columns:
        gset_dict[gset_name] = list(gset_df[gset_name].dropna())

    # Clean up (removes 1-4 genes per module)
    print('Absent from measured genes:\n',{key: [item for item in value if item not in adata.var_names] for key, value in gset_dict.items()})
    gset_dict = {key: np.array([item for item in value if item in adata.var_names]) for key, value in gset_dict.items()}

    gpres_dict = dict()
    output_file = f'Top200DEG_{leiRes}_{Experiment[9:-3]}.xlsx'
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        for gr in adata.obs[leiRes].cat.categories.to_list():
            #print(gr)
            df = sc.get.rank_genes_groups_df(adata, group=gr).head(200)

            df.to_excel(writer, sheet_name=gr, index=False)

            #print(df)
            gpres_dict[gr] = df[["names", "scores"]]

    dict_A = gset_dict.copy()
    dict_B = gpres_dict.copy()

    # Step 1: Create a DataFrame for heatmap
    all_genes = list(set(gene for genes in dict_A.values() for gene in genes))
    cell_types = sorted(list(dict_A.keys()))  # Sort cell types alphabetically

    # Initialize the DataFrame with NaN
    heatmap_data = pd.DataFrame(np.nan, index=all_genes, columns=cell_types)

    # Populate heatmap with scores from dict B for corresponding genes
    for cell_type, genes in dict_A.items():
        for gene in genes:
            heatmap_data.loc[gene, cell_type] = np.nan  # Initialize with NaN

    # Generate a heatmap for each group in dict B
    for group, df in dict_B.items():
        # Create a copy of the heatmap_data and populate with scores for the current group
        group_heatmap_data = heatmap_data.copy()
        
        for _, row in df.iterrows():
            gene = row['names']
            score = row['scores']
            for cell_type, genes in dict_A.items():
                if gene in genes:
                    group_heatmap_data.loc[gene, cell_type] = float(score)  # Cast score to float

        # Step 2: Fill NaN values with 0 (or another method)
        group_heatmap_data = group_heatmap_data.fillna(0)

        # Step 3: Filter heatmap data to only show genes with non-zero values
        filtered_heatmap_data = group_heatmap_data.loc[(group_heatmap_data != 0).any(axis=1)]

        # Step 4: Sort genes (Y-axis) by highest score, from high to low
        filtered_heatmap_data = filtered_heatmap_data.assign(max_score=filtered_heatmap_data.max(axis=1))
        filtered_heatmap_data = filtered_heatmap_data.sort_values(by='max_score', ascending=False)
        filtered_heatmap_data = filtered_heatmap_data.drop(columns='max_score')

        # Step 5: Adjust figure size dynamically based on the number of genes
        num_genes = len(filtered_heatmap_data)
        plt.figure(figsize=(10, max(5, num_genes * 0.5)))

        # Step 6: Standard heatmap with sorted X-axis and no grid lines
        ax = sns.heatmap(
            filtered_heatmap_data[cell_types],  # Ensure cell types are sorted on X-axis
            cmap="coolwarm",  # Color by scores
            cbar_kws={'label': 'z-score'},  # Colorbar label
            linewidths=0,  # Ensure grid is removed
            linecolor='none',
            square=False,
            yticklabels=True,  # Show gene labels
            xticklabels=True  # Show cell type labels
        )

        # Remove all spines and lines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        # Save the heatmap
        plt.title(f'{num_genes} from 200 Top DEG found in 2020 Signatures\nCluster {group} - {leiRes}')
        plt.savefig(f'figures/Heatmap_{group}_{leiRes}.png', bbox_inches='tight')  # Save with tight layout
        plt.close()



