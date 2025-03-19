import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
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

# Read new gne signature list
df = pd.read_csv(f'{refDir}/2024-08-20_gene_expression.csv')
# Pivot the DataFrame
new_gset_df = df.pivot(columns="Cell type ", values="Gene")
new_gset_df = new_gset_df.rename(columns={'Luminal 1': 'Epi_Luminal_1', 'Luminal 2': 'Epi_Luminal_2Psca'})

post_ext = '_post.h5ad'
lum_ext = '_annot_Lum_newSigs.h5ad'
all_ext = '_annot_All_newSigs.h5ad'

# ### Cell typing
# # From Wouter presentation
# # Cell types expected to be present: immune, mesenchymal, luminal, smooth muscle, basal (epi), neuroendocrine (epi), nerve
# # Mouse: ratio epi:stromal = 50:50
# # From Karthaus 2020 Science
# # Luminal 1/2/3 (all markers in presentation/paper), Glial, Myeloid, Mesenchymal 1/2, Vascular endothelium, Lymphatic endothelium, T/B/NK/dendritic (few), myofibroblasts/smooth muscle
# # Singulator protocol (data available right now) has different ratios of cell types. Much more L1 cells.
# # From Methods section:
# # splitting into epi (Epcam+), immune (Ptprc+), stromal (--). 
# ### big absences from 'measured_genes' (sample-specific): Cd19, Cd3*, Cd31/Pecam1, Sox10, some Wnt, Cd24a, ("Foxi1", "Atp6v1g3", "Atp6v1b1") discriminating the luminal_3, 

### Top Level Lineages (Epi,Imm,Strm) + SV removal
leiR = 3
leiRes = f'leiden_res{leiR}'
normToUse = 'SCT_data' #'scran_normalization', 'SCT_data', 'SCT_counts', 'log1p_norm', 'No_norm'

#adata import
adata = sc.read_h5ad(Experiment + post_ext)

# Remove SING bad samples: 'WK-1501_BL6_INTACT_AP_Test3_SORT' & 'WK-1384_BL6_Intact_AP_2_SLT'
adata = adata[~(np.isin(adata.obs['batch'], np.array(['WK-1501_BL6_INTACT_AP_Test3_SORT','WK-1384_BL6_Intact_AP_2_SLT'])))]

sc.tl.score_genes_cell_cycle(adata,
    s_genes=cc_df.iloc[:, 0].dropna().tolist(),
    g2m_genes=cc_df.iloc[:, 1].dropna().tolist())

sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_in_top_20_genes','pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'],multi_panel=True,groupby='batch',rotation=90,save=Experiment+'_PreQC.pdf')
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["batch"],save=Experiment+'_batch',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'_timePoint',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'_Counts_QC',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],save=Experiment+'_PctCounts_QC',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["scDblFinder_score","doublet_class"],save=Experiment+'_DbltClass_QC',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=['S_score','G2M_score','phase'],save=Experiment+'_CellCycle_QC',show=False)
sc.pl.umap(adata,color=leiRes,legend_loc="on data",save=f"{Experiment}_{leiRes}",title=Experiment,show=False)

#UMAP of Top level lineage scores

adata.obs['EpcamNorm'] = adata.layers[normToUse].T[adata.var.index == 'Epcam'][0]
adata.obs['PtprcNorm'] = adata.layers[normToUse].T[adata.var.index == 'Ptprc'][0]

sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='EpcamNorm',save=f'{Experiment}_Main_EpcamNormScore.png',title=f'{Experiment}_Epcam', cmap='RdBu_r',show=False,layer=normToUse, vmin=0, vmax=2)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='PtprcNorm',save=f'{Experiment}_Main_PtprcNormScore.png',title=f'{Experiment}_Ptprc', cmap='RdBu_r',show=False,layer=normToUse, vmin=0, vmax=2)

#True / False instead of score
adata.obs['EpcamExpressing'] = pd.Categorical(adata.obs['EpcamNorm'] > 0, categories=[True,False])
adata.obs['PtprcExpressing'] = pd.Categorical(adata.obs['PtprcNorm'] > 0, categories=[True,False])
#adata.obs['CD140aExpressing'] = pd.Categorical(adata.obs['CD140aNorm'] > 0, categories=[True,False])
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='EpcamExpressing',save=f'{Experiment}_Main_EpcamExpressing.png',title=f'{Experiment}_EpcamExpressing',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='PtprcExpressing',save=f'{Experiment}_Main_PtprcExpressing.png',title=f'{Experiment}_PtprcExpressing',show=False)

# Violin plots of Gene scores
# together
sc.pl.violin(adata,keys=['EpcamNorm','PtprcNorm'], groupby=leiRes,save=f'{Experiment}_Main_ScranScore_{leiRes}_violin.png',rotation=90,show=False)

### Top level lineage Annotation
HighPtprcClustersDict = get_highscore_clusters(adata,'PtprcNorm',leires=leiRes,medianthresh=0)
HighPtprcClusters = list(HighPtprcClustersDict.keys())

adata.obs['TopLevelLineages'] = 'NonImmune'
adata.obs.loc[[adata.obs[leiRes].isin(HighPtprcClusters)][0],'TopLevelLineages'] = 'Immune'
adata.obs['Annotation'] = pd.Categorical(['Other']*adata.n_obs,categories= ['Other'])
adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Immune'])
adata.obs.loc[[adata.obs[leiRes].isin(HighPtprcClusters)][0],'Annotation'] = 'Immune'

sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='TopLevelLineages',save=f'{Experiment}_Main_{leiRes}_TopLevelLineages.png',title=f'{Experiment}_{leiRes}_Top_Level_Lineages',show=False)
# Filter out Immune cells
adata_epi = adata[~adata.obs[leiRes].isin(HighPtprcClusters)].copy()
sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=leiRes,legend_loc="on data",save=f"{Experiment}_Main_Epi_{leiRes}.png",title=Experiment,show=False)


### Remove Basal
Basal_dict = {'Basal': list(new_gset_df['Basal'].dropna())}

#sc.tl.score_genes(adata_epi, ["Krt5", "Krt14", "Trp63"], ctrl_size=50, n_bins=25, score_name='Basal')
sc.tl.score_genes(adata_epi, list(Basal_dict.values())[0], ctrl_size=50, n_bins=25, score_name='Basal')

sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color='Basal',save=f'{Experiment}_Basal_EpiSubtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_epi,keys='Basal', groupby=leiRes,save=f'{Experiment}_Basal_ScoreGene_{leiRes}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

HighBasalClustersDict = get_highscore_clusters(adata_epi,'Basal',leires=leiRes,medianthresh=0.5)
HighBasalClusters = list(HighBasalClustersDict.keys())

adata_epi.obs['EpiLevelLineages'] = 'Luminal'
adata_epi.obs.loc[[adata_epi.obs[leiRes].isin(HighBasalClusters)][0],'EpiLevelLineages'] = 'Basal'
adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Basal'])
adata.obs.loc[[adata.obs[leiRes].isin(HighBasalClusters)][0],'Annotation'] = 'Basal'

sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color='EpiLevelLineages',save=f'{Experiment}_Main_{leiRes}_EpiLevelLineages.png',title=f'{Experiment}_{leiRes}_Epi_Level_Lineages',show=False)

#Subset adata_epi for Luminal clusters...
adata_lum = adata_epi[~adata_epi.obs[leiRes].isin(HighBasalClusters)]

# Redoing clustering in Luminal subset adata_epi: 
sc.pp.highly_variable_genes(adata_lum, layer=normToUse)
adata_lum.X = adata_lum.layers[normToUse]
sc.pp.pca(adata_lum, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_variance_ratio(adata_lum, log=True,save=Experiment+'_Epi'+'_PCA_Variance',show=False)
sc.pp.neighbors(adata_lum, use_rep="X_pca",n_neighbors=30, n_pcs=30)
sc.tl.umap(adata_lum)

leiR_Prec = 3
leiRes_Prec = f'leiden_res{leiR_Prec}'
sc.tl.leiden(adata_lum, key_added=leiRes_Prec, resolution=leiR_Prec)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=leiRes_Prec,legend_loc="on data",save=f"{Experiment}_Luminal_Recluster{leiRes_Prec}.png",title=Experiment,show=False)

### Remove SV
Lum_SV_dict = {'Epi_Luminal_SV': list(gset_df['Epi_Luminal_SV'].dropna())}
# Clean up 
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata_lum.var_names] for key, value in Lum_SV_dict.items()})
Lum_SV_dict = {key: np.array([item for item in value if item in adata_lum.var_names]) for key, value in Lum_SV_dict.items()}

sc.tl.score_genes(adata_lum, Lum_SV_dict['Epi_Luminal_SV'], ctrl_size=len(Lum_SV_dict['Epi_Luminal_SV']), n_bins=25, score_name='Luminal_SV')

sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color='Luminal_SV',save=f'{Experiment}_Luminal_SV_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_lum,keys='Luminal_SV', groupby=leiRes_Prec,save=f'{Experiment}_Luminal_SV_{leiRes_Prec}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

HighLSVClustersDict = get_highscore_clusters(adata_lum,'Luminal_SV',leires=leiRes_Prec)
M = list(HighLSVClustersDict.values())

HighLSVClusters = np.array(list(HighLSVClustersDict.keys()))[(np.median(M) + 6 * median_abs_deviation(M) < M)]

#Subset adata_lum for Luminal clusters...
adata_lumSV = adata_lum[adata_lum.obs[leiRes_Prec].isin(HighLSVClusters)].copy()
adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Luminal_SV'])
adata.obs.loc[adata.obs.index[np.isin(adata.obs.index,adata_lumSV.obs.index)],'Annotation'] = 'Luminal_SV'

adata_lum = adata_lum[~adata_lum.obs[leiRes_Prec].isin(HighLSVClusters)].copy()

### L1 / L2
luminal_marker_names = ['Epi_Luminal_1','Epi_Luminal_2Psca']
epi_gset_dict = dict()
for gset_name in luminal_marker_names:
    epi_gset_dict[gset_name] = list(new_gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata_lum.var_names] for key, value in epi_gset_dict.items()})
epi_gset_dict = {key: np.array([item for item in value if item in adata_lum.var_names]) for key, value in epi_gset_dict.items()}
epi_gset_dict['Epi_Luminal_1'] = np.append(epi_gset_dict['Epi_Luminal_1'], np.array('Gm42418'))

print('Removed common genes in L1&L2 signatures:')
print(np.intersect1d(sorted(epi_gset_dict['Epi_Luminal_1']),sorted(epi_gset_dict['Epi_Luminal_2Psca'])))

epi_gset_dict['Epi_Luminal_1'] = np.setdiff1d(epi_gset_dict['Epi_Luminal_1'], epi_gset_dict['Epi_Luminal_2Psca'])
epi_gset_dict['Epi_Luminal_2Psca'] = np.setdiff1d(epi_gset_dict['Epi_Luminal_2Psca'], epi_gset_dict['Epi_Luminal_1'])

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata_lum, gs, ctrl_size=len(gs), n_bins=25, score_name=gset_name,use_raw=False)
hist_gene_sigs(adata_lum,luminal_marker_names,'Luminal_Subtypes_ScoreGene',Experiment=Experiment,thresh_gs=None,log_bool=True)

allEpiscores = np.array([adata_lum.obs[gse].tolist() for gse in luminal_marker_names])
vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=luminal_marker_names,vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_Luminal_Subtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_lum,keys=luminal_marker_names, groupby=leiRes_Prec,save=f'{Experiment}_Luminal_ScoreGene_{leiRes_Prec}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

#Typing
diff_thresh = 0.10

# DIFF
typedeflist = []
df_work = adata_lum.obs[luminal_marker_names]
for row in df_work.iterrows():
    markers = ['L1','L2']

    scores = row[1:][0].tolist()
    maxInd = np.argmax(scores)

    maxType = markers.pop(maxInd)
    maxScore = scores.pop(maxInd)
    AmbigWith = [typ for i,typ in enumerate(markers) if maxScore - diff_thresh < scores[i]]
    
    if AmbigWith:
        typedeflist.append(f"Ambig_{'|'.join(sorted(AmbigWith+[maxType]))}")
    else:
        typedeflist.append(maxType)

adata_lum.obs[f'TypeDiffAmbig'] = pd.Categorical(typedeflist,categories=sorted(np.unique(typedeflist), key=len))

for c in adata_lum.obs[f'TypeDiffAmbig'].cat.categories:
    adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories([c])
    temp = adata_lum[adata_lum.obs[f'TypeDiffAmbig']==c].copy()
    adata.obs.loc[adata.obs.index[np.isin(adata.obs.index,temp.obs.index)],'Annotation'] = c

adata_lum.obs['TypeDiffAmbig'] = quant_singleSeries(adata_lum.obs['TypeDiffAmbig'])

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=f'TypeDiffAmbig',title=f'{Experiment} ScoreGene Lum Subtypes Min Diff Ambig+',save=f'{Experiment}_Luminal_Subtypes_ScoreGene_DiffAmbig.png',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=f'Annotation',title=f'{Experiment} ScoreGene Annotation',save=f'{Experiment}_Annotation.png',show=False)

sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["batch"],save=Experiment+'Lum_NoSV_batch',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'Lum_NoSV_timePoint',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'Lum_NoSV_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'Lum_NoSV_Counts_QC',show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],save=Experiment+'Lum_NoSV_PctCounts_QC',show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["scDblFinder_score","doublet_class"],save=Experiment+'Lum_NoSV_DbltClass_QC',show=False)

### BONUSES

# CEll cycle UMAPs
sc.tl.score_genes_cell_cycle(adata_lum,
    s_genes=cc_df.iloc[:, 0].dropna().tolist(),
    g2m_genes=cc_df.iloc[:, 1].dropna().tolist())
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=['S_score','G2M_score','phase'],save=Experiment+'Lum_NoSV_CellCycle_QC',show=False)

# Check All v2 scores
for typeG in new_gset_df.columns.tolist():
    typeG_genes = {typeG: list(new_gset_df[typeG].dropna())}
    # Clean up 
    print('Absent from measured genes:\n',{key: [item for item in value if item not in adata.var_names] for key, value in typeG_genes.items()})
    typeG_genes = {key: np.array([item for item in value if item in adata.var_names]) for key, value in typeG_genes.items()}
    sc.tl.score_genes(adata, typeG_genes[typeG], ctrl_size=len(typeG_genes[typeG]), n_bins=25, score_name=typeG+'_v2')
    sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=typeG+'_v2',save=f'{Experiment}_{typeG}_ScoreGene_v2.png', cmap='RdBu_r',show=False)


# Karthaus plotting 
# By L1 / L2
for tp in adata_lum.obs['timePoint'].cat.categories:
    karthaus_plotting(adata_lum,tp,Experiment,'TypeDiffAmbig','GeneScore')
# With Batch coloring
for tp in adata_lum.obs['timePoint'].cat.categories:
    adat = adata_lum[adata_lum.obs['timePoint']==tp].copy()
    if 'ENZ' in Experiment:
        adat = adat[adat.obs['tissueProv'].str.startswith('AP')]
    #adat.obs['batch'] = adat.obs['batch'].cat.add_categories(['Ambig'])
    #adat.obs.loc[adat.obs["TypeDiffAmbig"].str.startswith('Ambig'),'batch'] = 'Ambig'
    #adat.obs['test'] = adat.obs['batch'].copy()
    sc.pl.scatter(adat,size=400000/adat.n_obs,color='batch',x='Epi_Luminal_1',y='Epi_Luminal_2Psca',alpha=0.8,save=Experiment+'scattL1L2_Batch'+tp,show=False)

#How many cells are Luminal No SV
#HIST cells in samples, epi vs total 
plt.hist(adata.obs['batch'], bins=len(adata.obs['batch'].unique()), alpha=0.5, label='All')#, edgecolor='black')
plt.hist(adata_lum.obs['batch'], bins=len(adata_lum.obs['batch'].unique()), alpha=0.5, label='Luminal 1/2')#, edgecolor='black', linewidth=0.1)
plt.xlabel(None)
plt.ylabel('Number of cells')
plt.title(f'Luminal cells in {Experiment} Samples')
plt.legend()
plt.xticks(np.arange(0,len(adata_lum.obs['batch'].unique()),1)+.5,rotation=45, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.grid(False)
plt.savefig(f'./figures/{Experiment}_{leiRes}_ProportionLumNoSV_Barplot.pdf')
plt.close()

#Add annot column to adata_lum
adata_lum.obs['Annotation'] = adata.obs['Annotation'].reindex(adata_lum.obs.index)

# Writing Epi subset file
adata_lum.write(f'{Experiment}{lum_ext}')
adata.write(f'{Experiment}{all_ext}')


