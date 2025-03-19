##Annotate Cells

# Purpose:
#   Cell annotation / Cell typing
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import median_abs_deviation
import anndata as ad
import os
import warnings

from G03_Annotation_Functions import *

print(f'Scanpy: {sc.__version__}')
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'.*")

sc.set_figure_params(figsize=(6, 6), dpi_save=300)

### TO CHECK WHEN CHANGING SAMPLES ###
Experiment='Wouter21_SING'
### TO CHECK WHEN CHANGING SAMPLES ###

refDir = '/mnt/etemp/ahrmad/wouter/refs'
Dir = f'/mnt/etemp/ahrmad/wouter/{Experiment}'
if not (os.path.exists(refDir)) & (os.path.exists(Dir)):
    raise ValueError('Check folder paths')
gene_set_table = f'{refDir}/table_s8_summary.txt'
#Epi subtypes Supp gene lists
gset_df = pd.read_csv(gene_set_table, sep='\t')

post_ext = '_Post.h5ad'
#post_ext = '_GeneMat.h5ad'
epi_ext = '_annot_Epi.h5ad'
lum_ext = '_annot_Lum.h5ad'
all_ext = '_annot_All.h5ad'


#CREATION OF THE GENE MAT
#import snapatac2 as snap
#pm = snap.read(Experiment + '_Post.h5ad').to_memory()
## Create Gene matrix Anndata
## Build Gene matrix
#adata = snap.pp.make_gene_matrix(pm, snap.genome.mm10)
## Do some basic filtering 
#sc.pp.filter_genes(adata, min_cells= 5)
#sc.pp.normalize_total(adata)
#sc.pp.log1p(adata)
#sc.external.pp.magic(adata, solver="approximate")
#adata.obsm["X_umap"] = pm.obsm["X_umap"]

if 'ENZ' in Experiment:
    Epcam_median_threshold = 0.15
    Ptprc_median_threshold = 0.4
    Luminal_median_threshold_SV = 0.14
    Basal_median_threshold_SV = 0.11
    diff_thresh = 0.05
    min_score = 0
if 'SING' in Experiment:
    Epcam_median_threshold = 0.6
    Ptprc_median_threshold = 0.7
    Luminal_median_threshold_SV = 0.3
    Basal_median_threshold_SV = 0.25
    diff_thresh = 0.05
    min_score = 0

### Top Level Lineages (Epi,Imm,Strm) + SV removal
#adata import
if post_ext == '_Post.h5ad':
    import snapatac2 as snap
    adataSNAP = snap.read(Experiment + post_ext,backed=None)
    # Build Gene matrix
    gene_matrix = snap.pp.make_gene_matrix(adataSNAP, snap.genome.mm10)
    # Do some basic filtering 
    sc.pp.filter_genes(gene_matrix, min_cells= 5)
    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)
    adata = gene_matrix.copy()

    adata.obsm["X_umap"] = adataSNAP.obsm["X_umap"]
    del gene_matrix,adataSNAP

elif post_ext == '_GeneMat.h5ad':
    adata = sc.read_h5ad(Experiment + post_ext)

for qc_gene in ["Mt", "Rb", "Hb"]:
    adata.var[qc_gene] = adata.var.index.str.startswith(qc_gene)

sc.pp.calculate_qc_metrics(adata, qc_vars=["Mt", "Rb", "Hb"], inplace=True, percent_top=[20], log1p=True)
adata.obs['batch'] = pd.Categorical(adata.obs['batch'],categories=['WK-1501_BL6_INTACT_AP_Test3_SORT','WK-1384_BL6_Intact_AP_2_SLT', 'WK-1501_BL6_INTACT_AP_Test3','WK-1580_BL6_AP_RegenDay2', 'WK-1585_Castrate_Day28_AP_BL6','WK-1585_INTACT_AP_BL6_Citrate', 'WK-1585_INTACT_AP_BL6_Contrl','WK-1585_Regen_Day3_AP_BL6', 'WK-1580_BL6_AP_RegenDay1'])

sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_in_top_20_genes','pct_counts_Mt', 'pct_counts_Rb', 'pct_counts_Hb'],multi_panel=True,groupby='batch',rotation=90,save=Experiment+'_PreQC_MAGIC.pdf')

sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["batch"],save=Experiment+'_batch',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'_timePoint',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["n_fragment", "frac_dup", "frac_mito"],save=Experiment+'_Counts_QC',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["FRiP", "tsse"],save=Experiment+'_PctCounts_QC',show=False)
adata_ori = adata.copy()

#Redoing clustering with GENE associated: 
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_variance_ratio(adata, log=True,save=Experiment+'_PCA_Variance',show=False)
sc.pp.neighbors(adata, use_rep="X_pca",n_neighbors=30, n_pcs=20)
sc.tl.umap(adata)

sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["batch"],save=Experiment+'_batch',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'_timePoint',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["n_fragment", "frac_dup", "frac_mito"],save=Experiment+'_Counts_QC',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=["FRiP", "tsse"],save=Experiment+'_PctCounts_QC',show=False)


leiR = 5
leiRes = f'leiden_res{leiR}'

sc.tl.leiden(adata, key_added=leiRes, resolution=leiR)

sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=leiRes,legend_loc="on data",save=f"{Experiment}_Recluster{leiRes}.png",title=Experiment,show=False)

#sc.pl.umap(adata,color=leiRes,legend_loc="on data",save=f"{Experiment}_{leiRes}",title=Experiment,show=False)

#UMAP of Top level lineage scores
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=['Epcam','Ptprc'],save=f'{Experiment}_Main_Score.png', cmap='RdBu_r',show=False, vmin=0, vmax=2)
#True / False instead of score
#adata.obs['EpcamExpressing'] = pd.Categorical(adata.X.T[adata.var.index == 'Epcam'][0] > 0, categories=[True,False])
#sc.pl.umap(adata,color='EpcamExpressing',save=f'{Experiment}_Main_EpcamExpressing.png',title=f'{Experiment}_EpcamExpressing',size=300000/adata.n_obs,show=False)

# Violin plots of Gene scores
sc.pl.violin(adata,keys=['Epcam','Ptprc'], groupby=leiRes,save=f'{Experiment}_Main_ScranScore_{leiRes}_violin.png',rotation=90,show=False)

#adata.write(f'{Experiment}_GeneMatNOMAGIC.h5ad')
#adata = sc.read_h5ad(f'{Experiment}_GeneMatNOMAGIC.h5ad')


#Histograms of distribtuion of Top level lineage scores
#hist_SoloGenesScore(adata,["Epcam","Ptprc"],'Main_Epcam_Ptprc',Experiment=Experiment)
adata.obs['EpcamScore'] = adata.X.T[adata.var.index == 'Epcam'].toarray()[0]
adata.obs['PtprcScore'] = adata.X.T[adata.var.index == 'Ptprc'].toarray()[0]

#Get Epcam expressing clusters (> half the cells express Epcam)
HighEpcamClustersDict = get_highscore_clusters(adata,'EpcamScore',leires=leiRes,medianthresh=Epcam_median_threshold)
HighEpcamClusters = list(HighEpcamClustersDict.keys())
HighPtprcClustersDict = get_highscore_clusters(adata,'PtprcScore',leires=leiRes,medianthresh=Ptprc_median_threshold)
HighPtprcClusters = list(HighPtprcClustersDict.keys())
#Is there overlapping clusters with Immune
assert len(np.intersect1d(HighEpcamClusters,HighPtprcClusters))==0,f'Some clusters overlap between Epcam and Ptprc {np.intersect1d(HighEpcamClusters,HighPtprcClusters)}'

if 'ENZ' in Experiment:
    #Annotate adata
    adata.obs['TopLevelLineages'] = 'Stromal'
    adata.obs.loc[[adata.obs[leiRes].isin(HighPtprcClusters)][0],'TopLevelLineages'] = 'Immune'
    adata.obs.loc[[adata.obs[leiRes].isin(HighEpcamClusters)][0],'TopLevelLineages'] = 'Epithelial'
    sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='TopLevelLineages',save=f'{Experiment}_Main_{leiRes}_TopLevelLineages.png',title=f'{Experiment}_{leiRes}_Top_Level_Lineages',show=False)

    #Subset adata for Epcam clusters...
    adata_epi = adata[adata.obs[leiRes].isin(HighEpcamClusters)].copy()
    #adata_imm = adata[adata.obs[leiRes].isin(HighPtprcClusters)].copy()
    #adata_str = adata[~adata.obs[leiRes].isin(HighEpcamClusters+HighPtprcClusters)].copy()
    sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=leiRes,legend_loc="on data",save=f"{Experiment}_Main_Epi_{leiRes}.png",title=Experiment,show=False)

elif 'SING' in Experiment:
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
sc.tl.score_genes(adata_epi, ["Krt5", "Krt14", "Trp63"], ctrl_size=50, n_bins=25, score_name='Basal')

sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color='Basal',save=f'{Experiment}_Basal_EpiSubtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_epi,keys='Basal', groupby=leiRes,save=f'{Experiment}_Basal_ScoreGene_{leiRes}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

HighBasalClustersDict = get_highscore_clusters(adata_epi,'Basal',leires=leiRes,medianthresh=0)
HighBasalClusters = list(HighBasalClustersDict.keys())

adata_epi.obs['EpiLevelLineages'] = 'Luminal'
adata_epi.obs.loc[[adata_epi.obs[leiRes].isin(HighBasalClusters)][0],'EpiLevelLineages'] = 'Basal'

adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Basal'])
adata.obs.loc[[adata.obs[leiRes].isin(HighBasalClusters)][0],'Annotation'] = 'Basal'

sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color='EpiLevelLineages',save=f'{Experiment}_Main_{leiRes}_EpiLevelLineages.png',title=f'{Experiment}_{leiRes}_Epi_Level_Lineages',show=False)

#Subset adata_epi for Luminal clusters...
adata_lum = adata_epi[~adata_epi.obs[leiRes].isin(HighBasalClusters)].copy()

# Redoing clustering in Luminal subset adata_epi: 
sc.pp.highly_variable_genes(adata_lum)
sc.pp.pca(adata_lum, svd_solver="arpack", use_highly_variable=True)
#sc.pl.pca_variance_ratio(adata_lum, log=True,save=Experiment+'_Epi'+'_PCA_Variance',show=False)
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
adata_lum = adata_lum[~adata_lum.obs[leiRes_Prec].isin(HighLSVClusters)].copy()

adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Luminal_SV'])
adata.obs.loc[adata.obs.index[np.isin(adata.obs.index,adata_lumSV.obs.index)],'Annotation'] = 'Luminal_SV'

### L1L2
luminal_marker_names = ['Epi_Luminal_1','Epi_Luminal_2Psca']
epi_gset_dict = dict()
for gset_name in luminal_marker_names:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata_lum.var_names] for key, value in epi_gset_dict.items()})
epi_gset_dict = {key: np.array([item for item in value if item in adata_lum.var_names]) for key, value in epi_gset_dict.items()}

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
diff_thresh = 0.05

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
adata.write(f'{Experiment}{all_ext}')

#def assign_size(celltype):
#    return 220000/adata_lum.n_obs if celltype in ['Basal','L1','L2','L3'] else 120000/adata_lum.n_obs
#size_list = adata_lum.obs['TypeDiffAmbig'].apply(assign_size).tolist()


adata_lum.obs['TypeDiffAmbig'] = quant_singleSeries(adata_lum.obs['TypeDiffAmbig'])

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=f'TypeDiffAmbig',title=f'{Experiment} ScoreGene Lum Subtypes Min Diff Ambig+',save=f'{Experiment}_Luminal_Subtypes_ScoreGene_DiffAmbig.png',show=False)
#,size=size_list

#How many cells are Epi No SV
#HIST cells in samples, epi vs total 
if 'SING' in Experiment:
    plt.hist(adata.obs['batch'], bins=len(adata.obs['batch'].unique()), alpha=0.5, label='All')#, edgecolor='black')
    plt.hist(adata_lum.obs['batch'], bins=len(adata_lum.obs['batch'].unique()), alpha=0.5, label='Luminal 1/2')#, edgecolor='black', linewidth=0.1)
elif 'ENZ' in Experiment:
    plt.hist(adata.obs['batch'], bins=len(adata.obs['batch'].unique()), alpha=0.5, label='All', edgecolor='black')
    plt.hist(adata_lum.obs['batch'], bins=len(adata_lum.obs['batch'].unique()), alpha=0.5, label='Luminal 1/2', edgecolor='black', linewidth=0.3)  
plt.xlabel(None)
plt.ylabel('Number of cells')
plt.title(f'Luminal cells in {Experiment} Samples')
plt.legend()
plt.xticks(np.arange(0,len(adata_lum.obs['batch'].unique()),1)+.5,rotation=45, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.grid(False)
plt.savefig(f'{Experiment}_{leiRes}_ProportionLumNoSV_Barplot.pdf')
plt.close()

sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["batch"],save=Experiment+'Lum_NoSV_batch',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'Lum_NoSV_timePoint',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'Lum_NoSV_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'Lum_NoSV_Counts_QC',show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["n_fragment", "frac_dup", "frac_mito"],save=Experiment+'Lum_NoSV_Counts_QC',show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["FRiP", "tsse"],save=Experiment+'Lum_NoSV_FRiP_TSSe_QC',show=False)

for tp in adata_lum.obs['timePoint'].cat.categories:
    karthaus_plotting(adata_lum,tp,Experiment,'TypeDiffAmbig','GeneScore')

# Writing Epi subset file
adata_lum.write(f'{Experiment}{lum_ext}')



#Redoing clustering in Epi subset adata: 
sc.pp.highly_variable_genes(adata_epi)
sc.pp.pca(adata_epi, svd_solver="arpack", use_highly_variable=True)
#sc.pl.pca_variance_ratio(adata_epi, log=True,save=Experiment+'_Epi'+'_PCA_Variance',show=False)
sc.pp.neighbors(adata_epi, use_rep="X_pca",n_neighbors=30, n_pcs=20)
sc.tl.umap(adata_epi)

leiR_SV = 5
leiRes_SV = f'leiden_res{leiR_SV}'

sc.tl.leiden(adata_epi, key_added=leiRes_SV, resolution=leiR_SV)

sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=leiRes_SV,legend_loc="on data",save=f"{Experiment}_Main_Epi_Recluster{leiRes_SV}.png",title=Experiment,show=False)

epiSV_marker_names = ['Epi_Basal_SV','Epi_Luminal_SV']

epi_gset_dict = dict()
for gset_name in epiSV_marker_names:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
epi_gset_dict = {key: [item for item in value if item in adata_epi.var_names] for key, value in epi_gset_dict.items()}

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata_epi, gs, ctrl_size=len(gs), n_bins=25, score_name=gset_name,use_raw=False)

#hist_gene_sigs(adata_epi,epiSV_marker_names,'EpiSV_Subtypes_ScoreGene',Experiment=Experiment,thresh_gs=None,log_bool=True)

allEpiscores = np.array([adata_epi.obs[gse].tolist() for gse in epiSV_marker_names])
vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)
sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=epiSV_marker_names,vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_EpiSV_Subtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(16, 6), dpi_save=300)
sc.pl.violin(adata_epi,keys=epiSV_marker_names, groupby=leiRes_SV,save=f'{Experiment}_SV_ScoreGene_{leiRes_SV}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

# Same thing with ORA, the clusters agree....
#adata_epi_ora  = adata_epi.copy()
#adata_epi_ora.obs.drop(columns=['leiden_res1', 'Epi_Basal_SV', 'Epi_Luminal_SV'], inplace=True)
#
#markers_long_list = []
#for gset_name in epiSV_marker_names:
#    for gen in list(gset_df[gset_name].dropna()):
#        markers_long_list.append([gen,gset_name])
#
#markers_long = pd.DataFrame(markers_long_list)
#col_names = ['genesymbol', 'cell_type']
#markers_long.columns = col_names
#
#dc.run_ora(
#    mat=adata_epi_ora,
#    net=markers_long,
#    source='cell_type',
#    target='genesymbol',
#    min_n=10,
#    verbose=True,
#    use_raw=False
#)
#
#acts_epi = dc.get_acts(adata_epi_ora, obsm_key='ora_estimate')
## We need to remove inf and set them to the maximum value observed for pvals=0
#acts_epi_v = acts_epi.X.ravel()
#max_e = np.nanmax(acts_epi_v[np.isfinite(acts_epi_v)])
#acts_epi.X[~np.isfinite(acts_epi.X)] = max_e
#
#sc.pl.umap(acts_epi, color=epiSV_marker_names,save=f'{Experiment}_EpiSV_Subtypes_ORA.png',size=300000/acts_epi.n_obs, cmap='RdBu_r',show=False,vmin=0,vmax=20)
#
#sc.set_figure_params(figsize=(12, 6), dpi_save=300)
#sc.pl.violin(acts_epi,keys=epiSV_marker_names, groupby=leiRes_SV,save=f'{Experiment}_SV_ORA_{leiRes_SV}_violin.png',rotation=90,show=False)
#sc.set_figure_params(figsize=(6, 6), dpi_save=300)
#
#for subtype in epiSV_marker_names:
#    acts_epi.obs[subtype] = acts_epi.obsm['ora_estimate'][subtype].copy()
#hist_gene_sigs(acts_epi,epiSV_marker_names,'EpiSV_Subtypes_ORA',Experiment=Experiment,thresh_gs=None,log_bool=True)
#get_highscore_clusters(acts_epi,'Epi_Luminal_SV',leires=leiRes_SV,medianthresh=median_threshold_SV)
#get_highscore_clusters(acts_epi,'Epi_Basal_SV',leires=leiRes_SV,medianthresh=median_threshold_SV)

HighLuminalSVClusters = get_highscore_clusters(adata_epi,'Epi_Luminal_SV',leires=leiRes_SV,medianthresh=Luminal_median_threshold_SV)
HighBasalSVClusters = get_highscore_clusters(adata_epi,'Epi_Basal_SV',leires=leiRes_SV,medianthresh=Basal_median_threshold_SV)

adata_epi_NoSV = adata_epi[~adata_epi.obs[leiRes_SV].isin(HighLuminalSVClusters+HighBasalSVClusters)].copy()
adata_epi_NoSV = adata_epi

#How many cells are Epi No SV
#HIST cells in samples, epi vs total 
plt.hist(adata.obs['batch'], bins=len(adata.obs['batch'].unique()), alpha=0.5, label='All')
plt.hist(adata_epi_NoSV.obs['batch'], bins=len(adata_epi_NoSV.obs['batch'].unique()), alpha=0.5, label='Epithelial')
plt.xlabel(None)
plt.ylabel('Number of cells')
plt.title(f'Epithelial cells in {Experiment} Samples')
plt.legend()
plt.xticks(rotation=45, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.savefig(f'{Experiment}_Main_{leiRes}_ProportionEpiNoSV_Barplot.pdf')
plt.close()

sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["batch"],save=Experiment+'_EPI_NoSV_batch',title=Experiment,show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'_EPI_NoSV_timePoint',title=Experiment,show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'_EPI_NoSV_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["n_fragment", "frac_dup", "frac_mito"],save=Experiment+'_EPI_NoSV_Counts_QC',show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["FRiP", "tsse"],save=Experiment+'_EPI_NoSV_FRiP_TSSe_QC',show=False)

# Writing Epi subset file
adata_epi_NoSV.write(f'{Experiment}{epi_ext}')





### Finer Cell Typing (Basal, L1, L2 esp)

######################
######################
##score_genes METHOD##
######################
######################
##Method 1 with sc.tl.score_genes

#adata_epi_NoSV = sc.read_h5ad(f'{Experiment}{epi_ext}')

### Remove Basal
sc.tl.score_genes(adata_epi_NoSV, ["Krt5", "Krt14", "Trp63"], ctrl_size=50, n_bins=25, score_name='Basal')

sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color='Basal',save=f'{Experiment}_Basal_EpiSubtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_epi_NoSV,keys='Basal', groupby=leiRes,save=f'{Experiment}_Basal_ScoreGene_{leiRes}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

HighBasalClusters = get_highscore_clusters(adata_epi_NoSV,'Basal',leires=leiRes,medianthresh=0)

adata_epi_NoSV.obs['EpiLevelLineages'] = 'Luminal'
adata_epi_NoSV.obs.loc[[adata_epi_NoSV.obs[leiRes].isin(HighBasalClusters)][0],'EpiLevelLineages'] = 'Basal'

sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color='EpiLevelLineages',save=f'{Experiment}_Main_{leiRes}_EpiLevelLineages.png',title=f'{Experiment}_{leiRes}_Epi_Level_Lineages',show=False)

#Subset adata_epi for Luminal clusters...
adata_epi_NoSV = adata_epi_NoSV[~adata_epi_NoSV.obs[leiRes].isin(HighBasalClusters)].copy()
adata_epi_NoSV

epi_marker_names_NoSV = ['Epi_Luminal_1','Epi_Luminal_2Psca']

epi_gset_dict = dict()
for gset_name in epi_marker_names_NoSV:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata_epi_NoSV.var_names] for key, value in epi_gset_dict.items()})
epi_gset_dict = {key: np.array([item for item in value if item in adata_epi_NoSV.var_names]) for key, value in epi_gset_dict.items()}

print('Removed common genes in L1&L2 signatures:')
print(np.intersect1d(sorted(epi_gset_dict['Epi_Luminal_1']),sorted(epi_gset_dict['Epi_Luminal_2Psca'])))

epi_gset_dict['Epi_Luminal_1'] = np.setdiff1d(epi_gset_dict['Epi_Luminal_1'], epi_gset_dict['Epi_Luminal_2Psca'])
epi_gset_dict['Epi_Luminal_2Psca'] = np.setdiff1d(epi_gset_dict['Epi_Luminal_2Psca'], epi_gset_dict['Epi_Luminal_1'])

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata_epi_NoSV, gs, ctrl_size=len(gs), n_bins=25, score_name=gset_name,use_raw=False)

#hist_gene_sigs(adata_epi_NoSV,epi_marker_names_NoSV,'Epi_Subtypes_ScoreGene',Experiment=Experiment,thresh_gs=None,log_bool=False)

allEpiscores = np.array([adata_epi_NoSV.obs[gse].tolist() for gse in epi_marker_names_NoSV])
if 'ENZ' in Experiment:
    vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)
elif 'SING' in Experiment:
    vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)-1
    vminEpi,vmaxEpi = -0.2,0.6
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=epi_marker_names_NoSV,vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_Epi_NoSV_Subtypes_ScoreGene.png', cmap='RdBu_r',show=False)


#Typing

epi_marker_names_NoSV = ['Epi_Luminal_1','Epi_Luminal_2Psca']
# DIFF
typedeflist = []
df_work = adata_epi_NoSV.obs[epi_marker_names_NoSV]
for row in df_work.iterrows():
    markers = ['L1','L2']

    scores = row[1:][0].tolist()
    maxInd = np.argmax(scores)

    if maxInd < min_score:
        typedeflist.append('LowScore')
        continue

    maxType = markers.pop(maxInd)
    maxScore = scores.pop(maxInd)
    AmbigWith = [typ for i,typ in enumerate(markers) if maxScore - diff_thresh < scores[i]]
    
    if AmbigWith:
        typedeflist.append(f"Ambig_{'|'.join(sorted(AmbigWith+[maxType]))}")
    else:
        typedeflist.append(maxType)

adata_epi_NoSV.obs[f'TypeDiffAmbig'] = pd.Categorical(typedeflist,categories=sorted(np.unique(typedeflist), key=len))

adata_epi_NoSV.obs['TypeDiffAmbigNum'] = quant_singleSeries(adata_epi_NoSV.obs['TypeDiffAmbig'])

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=f'TypeDiffAmbigNum',title=f'{Experiment} ScoreGene Epi Subtypes Min Diff Ambig+',save=f'{Experiment}_Epi_Subtypes_ScoreGene_DiffAmbig{diff_thresh}.png',show=False)

for tp in adata_epi_NoSV.obs['timePoint'].cat.categories:
    karthaus_plotting(adata_epi_NoSV,tp,Experiment,'TypeDiffAmbig','GeneScore')

adata_epi_NoSV.write(f'{Experiment}{epi_ext[:-5]}_Subtypes_ScoreGene.h5ad')



















#QC
Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
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
dfs = []
for sample in Samples:
    atacdf = pd.read_csv(f'{refDir}/csv/{sample}_doublet_scores_ATACPeaks_Self.csv')
    gexdf = pd.read_csv(f'{refDir}/csv/{sample}_doublet_scores_CB_GEX.csv')
    mergDF = atacdf.merge(gexdf, on='obs_names', how='inner')
    
    mergDF['doublet_class'] = 'WHATAMI'
    mergDF.loc[(mergDF['doublet_class_x'] == 'singlet') & (mergDF['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
    mergDF.loc[(mergDF['doublet_class_x'] == 'doublet') & (mergDF['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet GEX Only'
    mergDF.loc[(mergDF['doublet_class_x'] == 'singlet') & (mergDF['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Singlet ATAC Only'
    mergDF.loc[(mergDF['doublet_class_x'] == 'doublet') & (mergDF['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
    
    assert ~any(mergDF['doublet_class'].str.contains('WHATAMI')),'Some cells have doublet status unassigned'
    assert mergDF['sample_x'].equals(mergDF['sample_y']),'Sample row not matching between GEX and ATAC'

    dfs.append(mergDF)

samples_df = pd.concat(dfs, ignore_index=True)
samples_df['cellnames'] = samples_df['obs_names'] + '_' + samples_df['sample_x']

samples_df_in_adata = samples_df[samples_df['cellnames'].isin(adata.obs.index)]

assert samples_df_in_adata['cellnames'].equals(adata.obs.index),'Cells in df and adata not in same order'

adata.obs['SCRUBLETdblt_class'] = pd.Categorical(samples_df_in_adata['doublet_class'])

sc.pl.umap(adata,color=["n_fragment", "tsse", "frac_dup", "frac_mito", "FRiP"],save=Experiment+'_UMAP_QC.pdf',show=False)
sc.pl.umap(adata,color=["SCRUBLETdblt_class"],palette=['blue','red'],save=Experiment+'_DbltSCRUBLET.pdf',show=False)



    

adata.var.index = adata.var.index.str.upper()

#Celltyping based on Supp gene lists
gset_df = pd.read_csv(gene_set_table, sep='\t')
gset_dict = dict()
for gset_name in gset_df.columns:
    gset_dict[gset_name] = gset_df[gset_name].dropna().str.upper().tolist()

gset_dict = {key: [item for item in value if item in adata.var_names] for key, value in gset_dict.items()}

# Cell Types UMAPs
for gset_name in gset_dict.keys():
    gs = gset_dict[gset_name]
    print(f'{gset_name}: {len(gs)}genes Karthaus2020 TableS8')
    
    ctrl = len(gs) if (len(gs) > 50) else 50 #genes
    sc.tl.score_genes(adata, gs, ctrl_size=ctrl, n_bins=25, score_name=gset_name,use_raw=False)

# Write 
adata.write(Experiment + annot_ext)

#Plotting
adata = sc.read_h5ad(Experiment + annot_ext)

#Luminal +basal gene scores distrib
marker_names = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Luminal_3Foxi1']

# Create a figure and three subplots
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
# Plot histograms for each column
for i, marker_name in enumerate(marker_names):
    axes[i].hist(adata.obs[marker_name], bins=100, color='skyblue')

    # Add labels and title
    axes[i].set_xlabel(marker_name + ' Score')
    axes[i].set_ylabel('Number of cells')
    axes[i].set_title(marker_name + ' Score Distribution')

plt.tight_layout()
plt.savefig(Experiment+'_ATACGeneScore_Hist.pdf')
plt.close()

# Cell Types UMAPs
for gset_name in gset_dict.keys():
    sc.pl.umap(adata,color=gset_name,save=f'{Experiment}_{gset_name}_AnnotSuppS8_ctrl{ctrl}.png',title=f'{Experiment}_{gset_name}',show=False)

# Cell Types Dotplot
def dotplotting(adata,gset_dic,leidenRes,Experimenti=Experiment):
    name = list(gset_dic.keys())[0][0:3]
    sc.pl.dotplot(
        adata,
        groupby=leidenRes,
        var_names=gset_dic,
        standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
        use_raw=False,
        title=f'{Experiment} {name}',
        save=f'{Experiment}_{name}_AnnotSuppS8_Dotplot.pdf',
        show=False
    )

gset_dict_dot_Imm = {key: [item for item in value] for (key, value) in gset_dict.items() if key.startswith('Imm')}
gset_dict_dot_Str = {key: [item for item in value] for (key, value) in gset_dict.items() if key.startswith('Str')}
gset_dict_dot_Epi = {key: [item for item in value] for (key, value) in gset_dict.items() if key.startswith('Epi')}

leidenRes = "leiden_res0.25"
dotplotting(adata,gset_dict_dot_Imm,leidenRes)
dotplotting(adata,gset_dict_dot_Str,leidenRes)
dotplotting(adata,gset_dict_dot_Epi,leidenRes)

# L1/L2 Signature
def karthaus_plotting(adata,timepoint,Experimenti=Experiment):
    x = adata.obs['Epi_Luminal_1']
    y = adata.obs['Epi_Luminal_2Psca']
    # Create a scatter plot with different colors for x and y values
    plt.scatter(x, y, s=8, marker='o',alpha=0.9)
    # Add labels and title
    plt.xlabel('L1 Signature Score')
    plt.ylabel('L2 Signature Score')
    plt.title(f'{Experimenti} {timepoint}')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"./figures/L1L2Sig_{Experimenti}_{timepoint}.pdf")
    plt.close()
    print(f"Plotted L1/L2 {Experimenti} {timepoint}")

if 'ENZ' in Experiment:
    i_adata = adata[adata.obs['batch'].str.contains('WK-1350_I')].copy()
    r3_adata = adata[adata.obs['batch'].str.contains('WK-1350_R3')].copy()

    karthaus_plotting(i_adata,'Intact')
    karthaus_plotting(r3_adata,'RegenDay3')

elif 'SING' in Experiment:

    i_adata = adata[adata.obs['batch'].str.upper().str.contains('INTACT')].copy()
    r3_adata = adata[adata.obs['batch'].str.contains('Day3')].copy()
    r2_adata = adata[adata.obs['batch'].str.contains('Day2')].copy()
    r1_adata = adata[adata.obs['batch'].str.contains('Day1')].copy()
    c28_adata = adata[adata.obs['batch'].str.contains('Day28')].copy()

    karthaus_plotting(i_adata,'Intact')
    karthaus_plotting(r3_adata,'RegenDay3')
    karthaus_plotting(r2_adata,'RegenDay2')
    karthaus_plotting(r1_adata,'RegenDay1')
    karthaus_plotting(c28_adata,'CastrateDay28')


#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Imm*.png --outfile ENZ_Imm_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Str*.png --outfile ENZ_Str_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Epi*.png --outfile ENZ_Epi_CB.pdf

#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Imm*.png --outfile SING_Imm_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Str*.png --outfile SING_Str_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Epi*.png --outfile SING_Epi_CB.pdf








































import scanpy as sc
import matplotlib.pyplot as plt
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
DirExp = '/mnt/etemp/ahrmad/wouter/ENZ_ATAC'
qc_ext = '_filt.h5ad'
refDir = '/mnt/etemp/ahrmad/wouter/refs'


# Annotation
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






