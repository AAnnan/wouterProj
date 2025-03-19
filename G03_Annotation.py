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

# Read the data into a DataFrame
df = pd.read_csv(f'{refDir}/2024-08-20_gene_expression.csv')
# Pivot the DataFrame
new_gset_df = df.pivot(columns="Cell type ", values="Gene")
new_gset_df = new_gset_df.rename(columns={'Luminal 1': 'Epi_Luminal_1', 'Luminal 2': 'Epi_Luminal_2Psca'})

post_ext = '_post.h5ad'
lum_ext = '_annot_Lum_20240821.h5ad'
all_ext = '_annot_All.h5ad'

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

if 'ENZ' in Experiment:
    leiR = 4
    leiRes = f'leiden_res{leiR}'
if 'SING' in Experiment:
    leiR = 3
    leiRes = f'leiden_res{leiR}'

normToUse = 'SCT_data' #'scran_normalization', 'SCT_data', 'SCT_counts', 'log1p_norm', 'No_norm'

#adata import
adata = sc.read_h5ad(Experiment + post_ext)

# Remove SING bad samples: 'WK-1501_BL6_INTACT_AP_Test3_SORT' & 'WK-1384_BL6_Intact_AP_2_SLT'
if 'SING' in Experiment:
    adata = adata[~(np.isin(adata.obs['batch'], np.array(['WK-1501_BL6_INTACT_AP_Test3_SORT','WK-1384_BL6_Intact_AP_2_SLT'])))]
# Remove ENZ non-AP samples:
if 'ENZ' in Experiment:
    adata = adata[adata.obs['tissueProv'].str.startswith('AP')]
    sc.tl.leiden(adata, key_added=leiRes, resolution=leiR)

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


mouse_stromal_cell_markers = {
    'Vim': 'Vimentin',
    'Des': 'Desmin',
    'Pdgfra': 'Platelet-Derived Growth Factor Receptor Alpha',
    'Thy1': 'CD90',
    'Fn1': 'Fibronectin',
    'Col1a1': 'Collagen Type I Alpha 1 Chain'}
#
#mouse_stromal_cell_markers = {
#    'Vim': 'Vimentin',
#    'Fap': 'Fibroblast Activation Protein',
#    'Acta2': 'Alpha-Smooth Muscle Actin',
#    'Des': 'Desmin',
#    'Pdgfra': 'Platelet-Derived Growth Factor Receptor Alpha',
#    'Thy1': 'CD90',
#    'Nt5e': 'CD73',
#    'Eng': 'CD105',
#    'Fn1': 'Fibronectin',
#    'Hapln1': 'Hyaluronan and Proteoglycan Link Protein 1',
#    'Postn': 'Periostin',
#    'Col1a1': 'Collagen Type I Alpha 1 Chain'}
#
#for g,d in mouse_stromal_cell_markers.items():
#    print(f'{g}: {d}')
#    print(f'{np.sum(adata.layers[normToUse].T[adata.var.index == g][0] > 0)}/{adata.n_obs}')
#
#mouse_epi_cell_markers = {
#    'Krt5': 'keratin-5 ',
#    'Krt8': 'keratin-8 ',
#    'Krt18': 'keratin-18',
#    'Cdh1': 'E-Cadherin '}
#
#for g,d in mouse_epi_cell_markers.items():
#    print(f'{g}: {d}')
#    print(f'{np.sum(adata.layers[normToUse].T[adata.var.index == g][0] > 0)}/{adata.n_obs}')
#
#sc.tl.score_genes(adata, mouse_epi_cell_markers.keys(), ctrl_size=50, n_bins=25, score_name='EpiScore')
#
##Histograms of distribtuion of Top level lineage scores
#
#hist_SoloGenesScore(adata,["Epcam","Ptprc"],'Main_Epcam_Ptprc',Experiment=Experiment)
#hist_SoloGenesScore(adata,["Epcam","Ptprc"],'Main_Epcam_Ptprc_raw',Experiment=Experiment,raw=True)
#hist_SoloGenesScore(adata,["EpiScore"],'Main_Episcore',Experiment=Experiment)
#hist_SoloGenesScore(adata,mouse_epi_cell_markers.keys(),'Main_EpiScore_raw',Experiment=Experiment,raw=True)
#hist_SoloGenesScore(adata,mouse_stromal_cell_markers.keys(),'Main_StrScore_raw',Experiment=Experiment,raw=True)
#hist_SoloGenesScore(adata,mouse_stromal_cell_markers1.keys(),'Main_StrScore1_raw',Experiment=Experiment,raw=True)
#hist_SoloGenesScore(adata,['Prom1','Pbsn','Ly6a','Nkx3-1'],'OtherGenes_raw',Experiment=Experiment,raw=True)
#
#plt.hist(adata.obs['EpiScore'],bins=100)
#plt.savefig(Experiment+'_EpiScore.png')
#plt.close()
#
#gene_pres_Wouter = ['Epcam','Krt5','Krt8','Krt18','Prom1','Pbsn','Svs2','Svs5','Ly6a','Nkx3-1']
#testing = False
#if testing:
#    adata_2 = adata[np.isin(adata.obs['batch'], np.array(['WK-1585_INTACT_AP_BL6_Citrate','WK-1585_INTACT_AP_BL6_Contrl']))].copy()
#    sc.pp.highly_variable_genes(adata_2, layer=normToUse)
#    adata_2.X = adata_2.layers[normToUse]
#    sc.pp.pca(adata_2, svd_solver="arpack", use_highly_variable=True)
#    sc.pp.neighbors(adata_2, use_rep="X_pca",n_neighbors=20, n_pcs=20)
#    sc.tl.umap(adata_2)
#    sc.pl.umap(adata_2,alpha=0.8,color=gene_pres_Wouter,color_map='jet',layer=normToUse,legend_loc="on data",save=f"{Experiment}_gene_pres_Wouter.png",show=False)
#    sc.pl.umap(adata_2,size=400000/adata_2.n_obs,alpha=0.8,color='batch',save=f"{Experiment}_gene_pres_Wouter_batch.png",show=False)

if normToUse=='SCT_data':
    adata.obs['EpcamNorm'] = adata.layers[normToUse].T[adata.var.index == 'Epcam'][0]
    adata.obs['PtprcNorm'] = adata.layers[normToUse].T[adata.var.index == 'Ptprc'][0]
    #adata.obs['CD140aNorm'] = adata.layers[normToUse].T[adata.var.index == 'Pdgfra'][0]
elif normToUse=='SCT_counts':
    adata.X = adata.layers[normToUse].copy()
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers[normToUse+"log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    normToUse = normToUse+"log1p_norm"
    adata.obs['EpcamNorm'] = adata.layers[normToUse].T[adata.var.index == 'Epcam'][0]
    adata.obs['PtprcNorm'] = adata.layers[normToUse].T[adata.var.index == 'Ptprc'][0]
else:
    adata.obs['EpcamNorm'] = adata.layers[normToUse].T[adata.var.index == 'Epcam'].toarray()[0]
    adata.obs['PtprcNorm'] = adata.layers[normToUse].T[adata.var.index == 'Ptprc'].toarray()[0]

#UMAP of Top level lineage scores
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='EpcamNorm',save=f'{Experiment}_Main_EpcamNormScore.png',title=f'{Experiment}_Epcam', cmap='RdBu_r',show=False,layer=normToUse, vmin=0, vmax=2)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='PtprcNorm',save=f'{Experiment}_Main_PtprcNormScore.png',title=f'{Experiment}_Ptprc', cmap='RdBu_r',show=False,layer=normToUse, vmin=0, vmax=2)
#sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='StromalScore',save=f'{Experiment}_Main_StromalScore.png',title=f'{Experiment}_StromalScore', cmap='RdBu_r',show=False, vmin=0, vmax=0.5)
#sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='EpiScore',save=f'{Experiment}_Main_EpiScore.png',title=f'{Experiment}_EpiScore', cmap='RdBu_r',show=False, vmin=0, vmax=0.5)
#sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='CD140aNorm',save=f'{Experiment}_Main_CD140aNormScore.png',title=f'{Experiment}_CD140a', cmap='RdBu_r',show=False,layer=normToUse, vmin=0, vmax=2)
#True / False instead of score
adata.obs['EpcamExpressing'] = pd.Categorical(adata.obs['EpcamNorm'] > 0, categories=[True,False])
adata.obs['PtprcExpressing'] = pd.Categorical(adata.obs['PtprcNorm'] > 0, categories=[True,False])
#adata.obs['CD140aExpressing'] = pd.Categorical(adata.obs['CD140aNorm'] > 0, categories=[True,False])
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='EpcamExpressing',save=f'{Experiment}_Main_EpcamExpressing.png',title=f'{Experiment}_EpcamExpressing',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='PtprcExpressing',save=f'{Experiment}_Main_PtprcExpressing.png',title=f'{Experiment}_PtprcExpressing',show=False)
#sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='CD140aExpressing',save=f'{Experiment}_Main_CD140aExpressing.png',title=f'{Experiment}_CD140aExpressing',show=False)
#True / False Stromal [PTPRC-,ETPCAM-,CD140a+]
#adata.obs['StromalExpressing'] = pd.Categorical((adata.obs['EpcamNorm'] <= 0) & (adata.obs['PtprcNorm'] <= 0) & (adata.obs['CD140aNorm'] > 0), categories=[True,False])
#sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='StromalExpressing',save=f'{Experiment}_Main_StromalExpressing.png',title=f'{Experiment} Stromal\n[PTPRC-,ETPCAM-,CD140a+]',show=False)
#adata.obs['StromalExpressing2'] = pd.Categorical((adata.obs['EpcamNorm'] <= 0) & (adata.obs['PtprcNorm'] <= 0), categories=[True,False])
#sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='StromalExpressing2',save=f'{Experiment}_Main_StromalExpressing2.png',title=f'{Experiment} Stromal\n[PTPRC-,ETPCAM-]',show=False)

# Violin plots of Gene scores
# To sort
#mean_values_Ep = adata.obs.groupby(leiRes,observed=False)['EpcamNorm'].median()
#sorted_clusters_Ep = mean_values_Ep.sort_values(ascending=False).index.tolist()
#mean_values_Pt = adata.obs.groupby(leiRes,observed=False)['PtprcNorm'].median()
#sorted_clusters_Pt = mean_values_Pt.sort_values(ascending=False).index.tolist()
# add this to the violin plot parameters: order=sorted_clusters_Pt

#sc.pl.violin(adata,keys='PtprcNorm', groupby=leiRes,save=f'{Experiment}_Main_PtprcNormScore_{leiRes}_violin.png',order=sorted_clusters_Pt,show=False)
# together
sc.pl.violin(adata,keys=['EpcamNorm','PtprcNorm'], groupby=leiRes,save=f'{Experiment}_Main_ScranScore_{leiRes}_violin.png',rotation=90,show=False)
# old box plots
#boxplot_gene_sigs(adata,'EpcamNorm','RNA',clusters=leiRes)
#boxplot_gene_sigs(adata,'PtprcNorm','RNA',clusters=leiRes)

if 'ENZ' in Experiment:
    sc.tl.score_genes(adata, mouse_stromal_cell_markers.keys(), ctrl_size=50, n_bins=25, score_name='StromalScore')
    sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='StromalScore',save=f'{Experiment}_Main_StromalExpressing.png',title=f'{Experiment}_StromalExpressing',show=False)
    
    HighStromalClustersDict = get_highscore_clusters(adata,'StromalScore',leires=leiRes,medianthresh=0)
    HighStromalClusters = list(HighStromalClustersDict.keys())

    #Get Epcam expressing clusters (> half the cells express Epcam)
    HighEpcamClustersDict = get_highscore_clusters(adata,'EpcamNorm',leires=leiRes,medianthresh=0)
    HighEpcamClusters = list(HighEpcamClustersDict.keys())

    HighPtprcClustersDict = get_highscore_clusters(adata,'PtprcNorm',leires=leiRes,medianthresh=0)
    HighPtprcClusters = list(HighPtprcClustersDict.keys())

    #Is there overlapping clusters with Immune
    if len(np.intersect1d(HighEpcamClusters,HighPtprcClusters,HighStromalClusters))>0:
        print('Some clusters overlap between Epcam, Ptprc and Stromal:')
        print(np.intersect1d(HighEpcamClusters,HighPtprcClusters,HighStromalClusters))
        print('Choosing highest median cluster type')
        for clu in list((set(HighEpcamClusters) & set(HighPtprcClusters)) | (set(HighPtprcClusters) & set(HighStromalClusters)) | (set(HighStromalClusters) & set(HighEpcamClusters))):
            for my_dict in [HighPtprcClustersDict,HighStromalClustersDict,HighEpcamClustersDict]:
                if clu not in my_dict:
                    my_dict[clu] = -100
            if HighEpcamClustersDict[clu] == max(HighPtprcClustersDict[clu],HighStromalClustersDict[clu],HighEpcamClustersDict[clu]):
                HighPtprcClustersDict.pop(clu)
                HighStromalClustersDict.pop(clu)
                print("EPCAM")
            elif HighPtprcClustersDict[clu] == max(HighPtprcClustersDict[clu],HighStromalClustersDict[clu],HighEpcamClustersDict[clu]):
                HighEpcamClustersDict.pop(clu)
                HighStromalClustersDict.pop(clu)
                print("PTPRC")
            elif HighStromalClustersDict[clu] == max(HighPtprcClustersDict[clu],HighStromalClustersDict[clu],HighEpcamClustersDict[clu]):
                HighEpcamClustersDict.pop(clu)
                HighPtprcClustersDict.pop(clu)
                print("Stromal")
            else:
                exit('Some clusters STILL overlap between Epcam, Ptprc and Stromal')
            
        HighEpcamClusters = list(HighEpcamClustersDict.keys())
        HighPtprcClusters = list(HighPtprcClustersDict.keys())
        HighStromalClusters = list(HighStromalClustersDict.keys())

    #Annotate adata
    adata.obs['TopLevelLineages'] = 'Other_Ambig'
    adata.obs.loc[[adata.obs[leiRes].isin(HighPtprcClusters)][0],'TopLevelLineages'] = 'Immune'
    adata.obs.loc[[adata.obs[leiRes].isin(HighEpcamClusters)][0],'TopLevelLineages'] = 'Epithelial'
    adata.obs.loc[[adata.obs[leiRes].isin(HighStromalClusters)][0],'TopLevelLineages'] = 'Stromal'
    sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='TopLevelLineages',save=f'{Experiment}_Main_{leiRes}_TopLevelLineages.png',title=f'{Experiment}_{leiRes}_Top_Level_Lineages',show=False)

    adata.obs['Annotation'] = pd.Categorical(['Other_Ambig']*adata.n_obs,categories= ['Other_Ambig'])
    adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Immune'])
    adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Epithelial'])
    adata.obs['Annotation'] = adata.obs['Annotation'].cat.add_categories(['Stromal'])
    adata.obs.loc[[adata.obs[leiRes].isin(HighPtprcClusters)][0],'Annotation'] = 'Immune'
    adata.obs.loc[[adata.obs[leiRes].isin(HighEpcamClusters)][0],'Annotation'] = 'Epithelial'
    adata.obs.loc[[adata.obs[leiRes].isin(HighStromalClusters)][0],'Annotation'] = 'Stromal'
    

    #Subset adata for Epcam clusters...
    adata_epi = adata[adata.obs[leiRes].isin(HighEpcamClusters)].copy()
    #adata_imm = adata[adata.obs[leiRes].isin(HighPtprcClusters)].copy()
    #adata_str = adata[~adata.obs[leiRes].isin(HighEpcamClusters+HighPtprcClusters)].copy()
    sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=leiRes,legend_loc="on data",save=f"{Experiment}_Main_Epi_{leiRes}.png",title=Experiment,show=False)

elif 'SING' in Experiment:
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
adata_lum = adata_epi[~adata_epi.obs[leiRes].isin(HighBasalClusters)].copy()

# Redoing clustering in Luminal subset adata_epi: 
sc.pp.highly_variable_genes(adata_lum, layer=normToUse)
adata_lum.X = adata_lum.layers[normToUse]
sc.pp.pca(adata_lum, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_variance_ratio(adata_lum, log=True,save=Experiment+'_Epi'+'_PCA_Variance',show=False)

if 'ENZ' in Experiment:
    sc.pp.neighbors(adata_lum, use_rep="X_pca",n_neighbors=30, n_pcs=10)
elif 'SING' in Experiment:
    sc.pp.neighbors(adata_lum, use_rep="X_pca",n_neighbors=30, n_pcs=30)
sc.tl.umap(adata_lum)
leiR_Prec = 3
leiRes_Prec = f'leiden_res{leiR_Prec}'
sc.tl.leiden(adata_lum, key_added=leiRes_Prec, resolution=leiR_Prec)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=leiRes_Prec,legend_loc="on data",save=f"{Experiment}_Luminal_Recluster{leiRes_Prec}.png",title=Experiment,show=False)

if 'SING' in Experiment:
    Stromal = {'Stromal': list(new_gset_df['Stromal'].dropna())}
    # Clean up 
    print('Absent from measured genes:\n',{key: [item for item in value if item not in adata.var_names] for key, value in Stromal.items()})
    Stromal = {key: np.array([item for item in value if item in adata.var_names]) for key, value in Stromal.items()}

    sc.tl.score_genes(adata, Stromal['Stromal'], ctrl_size=len(Stromal['Stromal']), n_bins=25, score_name='Stromal')

    sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color='Stromal',save=f'{Experiment}_Stromal_ScoreGene_v2.png', cmap='RdBu_r',show=False)

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

#Check Stromal scores
#sc.tl.score_genes(adata_lum, mouse_stromal_cell_markers.keys(), ctrl_size=50, n_bins=25, score_name='StromalScore')
#sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color='StromalScore',save=f'{Experiment}_Lum_StromalScore.png',title=f'{Experiment}_StromalScore', cmap='RdBu_r',show=False)


### L1L2


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

#def assign_size(celltype):
#    return 220000/adata_lum.n_obs if celltype in ['Basal','L1','L2','L3'] else 120000/adata_lum.n_obs
#size_list = adata_lum.obs['TypeDiffAmbig'].apply(assign_size).tolist()


adata_lum.obs['TypeDiffAmbig'] = quant_singleSeries(adata_lum.obs['TypeDiffAmbig'])

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=f'TypeDiffAmbig',title=f'{Experiment} ScoreGene Lum Subtypes Min Diff Ambig+',save=f'{Experiment}_Luminal_Subtypes_ScoreGene_DiffAmbig.png',show=False)
sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=f'Annotation',title=f'{Experiment} ScoreGene Annotation',save=f'{Experiment}_Annotation.png',show=False)
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
plt.savefig(f'./figures/{Experiment}_{leiRes}_ProportionLumNoSV_Barplot.pdf')
plt.close()

sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["batch"],save=Experiment+'Lum_NoSV_batch',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'Lum_NoSV_timePoint',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'Lum_NoSV_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'Lum_NoSV_Counts_QC',show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],save=Experiment+'Lum_NoSV_PctCounts_QC',show=False)
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=["scDblFinder_score","doublet_class"],save=Experiment+'Lum_NoSV_DbltClass_QC',show=False)

sc.tl.score_genes_cell_cycle(adata_lum,
    s_genes=cc_df.iloc[:, 0].dropna().tolist(),
    g2m_genes=cc_df.iloc[:, 1].dropna().tolist())
sc.pl.umap(adata_lum,size=400000/adata_lum.n_obs,alpha=0.8,color=['S_score','G2M_score','phase'],save=Experiment+'Lum_NoSV_CellCycle_QC',show=False)


# KArthaus plotting 
# By L1 / L2
for tp in adata_lum.obs['timePoint'].cat.categories:
    karthaus_plotting(adata_lum,tp,Experiment,'TypeDiffAmbig','GeneScore')

# By Batch
for tp in adata_lum.obs['timePoint'].cat.categories:
    adat = adata_lum[adata_lum.obs['timePoint']==tp].copy()

    if 'ENZ' in Experiment:
        adat = adat[adat.obs['tissueProv'].str.startswith('AP')]

    #adat.obs['batch'] = adat.obs['batch'].cat.add_categories(['Ambig'])
    #adat.obs.loc[adat.obs["TypeDiffAmbig"].str.startswith('Ambig'),'batch'] = 'Ambig'
    #adat.obs['test'] = adat.obs['batch'].copy()
    sc.pl.scatter(adat,size=400000/adat.n_obs,color='batch',x='Epi_Luminal_1',y='Epi_Luminal_2Psca',alpha=0.8,save=Experiment+'scattL1L2_Batch'+tp,show=False)

# Writing Epi subset file
adata_lum.write(f'{Experiment}{lum_ext}')

adata.write(f'{Experiment}{all_ext}')







































normToUse = 'No_norm'
adata_lum = sc.read_h5ad(f'{Experiment}{lum_ext}')
#l1l2 = [gene for sublist in epi_gset_dict.values() for gene in sublist]
adata_lum_GeneCheckL1 = adata_lum[:, epi_gset_dict['Epi_Luminal_1']].copy()
adata_lum_GeneCheckL2 = adata_lum[:, epi_gset_dict['Epi_Luminal_2Psca']].copy()

sc.set_figure_params(figsize=(17, 4), dpi_save=300)
#sc.pl.violin(adata_lum_GeneCheckL1,keys=epi_gset_dict['Epi_Luminal_1'],save=f'{Experiment}_L1GeneSet_violin.png',rotation=90,show=False)
#sc.pl.violin(adata_lum_GeneCheckL2,keys=epi_gset_dict['Epi_Luminal_2Psca'],save=f'{Experiment}_L2GeneSet_violin.png',rotation=90,show=False)
sc.pl.violin(adata_lum_GeneCheckL1,keys=epi_gset_dict['Epi_Luminal_1'],layer=normToUse,save=f'{Experiment}_L1GeneSet{normToUse}_violin.png',rotation=90,show=False)
sc.pl.violin(adata_lum_GeneCheckL2,keys=epi_gset_dict['Epi_Luminal_2Psca'],layer=normToUse,save=f'{Experiment}_L2GeneSet{normToUse}_violin.png',rotation=90,show=False)



# Stromal typing
str_df = gset_df.filter(like='Str_', axis=1)

str_gset_dict = dict()
for gset_name in str_df.columns:
    str_gset_dict[gset_name] = list(str_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata.var_names] for key, value in str_gset_dict.items()})
str_gset_dict = {key: np.array([item for item in value if item in adata.var_names]) for key, value in str_gset_dict.items()}

# Gene scoring on epi subset
for gset_name in str_gset_dict.keys():
    gs = str_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata, gs, ctrl_size=len(gs), n_bins=25, score_name=gset_name,use_raw=False)

sc.pl.umap(adata,size=400000/adata.n_obs,alpha=0.8,color=str_gset_dict.keys(),save=f'{Experiment}_StromalSignatures.png', cmap='RdBu_r',show=False, vmin=0, vmax=0.5)















######################
######################
######ORA METHOD######
######################
######################
#Method 2 with dc.run_ora

if 'SING' in Experiment:
    adata_lum = adata_lum[~(adata_lum.obs['batch']=='WK-1501_BL6_INTACT_AP_Test3_SORT')].copy()

import decoupler as dc

markers_long_list = []
for gset_name in epi_gset_dict:
    for gen in list(epi_gset_dict[gset_name]):
        markers_long_list.append([gen,gset_name])

markers_long = pd.DataFrame(markers_long_list,columns = ['genesymbol', 'cell_type'])
epi_marker_names_NoSV = ['Epi_Luminal_1','Epi_Luminal_2Psca']

dc.run_ora(
    mat=adata_lum,
    net=markers_long,
    source='cell_type',
    target='genesymbol',
    min_n=10,
    verbose=True,
    use_raw=False
)

acts_epi_NoSV = dc.get_acts(adata_lum, obsm_key='ora_estimate')
# We need to remove inf and set them to the maximum value observed for pvals=0
acts_epi_v = acts_epi_NoSV.X.ravel()
max_e = np.nanmax(acts_epi_v[np.isfinite(acts_epi_v)])
acts_epi_NoSV.X[~np.isfinite(acts_epi_NoSV.X)] = max_e

df_ora = acts_epi_NoSV.obsm['ora_estimate']
df_genescore = adata_lum.obs[epi_marker_names_NoSV]

assert df_ora.index.equals(df_genescore.index)

# Create a scatter plot with different colors for x and y values
plt.scatter(df_ora['Epi_Luminal_1'], df_genescore['Epi_Luminal_1'], s=8,label=f'Pearson cor = {round(df_ora['Epi_Luminal_1'].corr(df_genescore['Epi_Luminal_1']),2)}')

# Add labels and title
plt.xlabel('L1 Signature Score ORA')
plt.ylabel('L1 Signature Score AddModuleScore')
plt.title(f'{Experiment}_L1_ORA_AddModuleScore')
plt.legend()
plt.grid(False)
plt.tight_layout()
plt.savefig(f"{Experiment}_L1_ORA_AddModuleScore.pdf")
plt.close()

# Create a scatter plot with different colors for x and y values
plt.scatter(df_ora['Epi_Luminal_2Psca'], df_genescore['Epi_Luminal_2Psca'], s=8,label=f'Pearson cor = {round(df_ora['Epi_Luminal_2Psca'].corr(df_genescore['Epi_Luminal_2Psca']),2)}')

# Add labels and title
plt.xlabel('L2 Signature Score ORA')
plt.ylabel('L2 Signature Score AddModuleScore')
plt.title(f'{Experiment}_L2_ORA_AddModuleScore')
plt.legend()
plt.grid(False)
plt.tight_layout()
plt.savefig(f"{Experiment}_L2_ORA_AddModuleScore.pdf")
plt.close()


#Typing
diff_thresh = 3
step_thresh,max_thresh,min_thresh=1,50,4

# DIFF
typedeflist = []
df_work = acts_epi_NoSV.obsm['ora_estimate'][epi_marker_names_NoSV]

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

acts_epi_NoSV.obs[f'TypeDiffAmbig'] = pd.Categorical(typedeflist,categories=sorted(np.unique(typedeflist), key=len))

acts_epi_NoSV.obs['TypeDiffAmbig'] = quant_singleSeries(acts_epi_NoSV.obs['TypeDiffAmbig'])

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(acts_epi_NoSV,size=400000/acts_epi_NoSV.n_obs,alpha=0.8,color=f'TypeDiffAmbig',title=f'{Experiment} ORA Epi Subtypes Min Diff Ambig+',save=f'{Experiment}_Epi_Subtypes_ORA_DiffAmbig.png',show=False)
#,size=size_list

acts_epi_NoSV.obs['Epi_Luminal_1'] = acts_epi_NoSV.obsm['ora_estimate']['Epi_Luminal_1']
acts_epi_NoSV.obs['Epi_Luminal_2Psca'] = acts_epi_NoSV.obsm['ora_estimate']['Epi_Luminal_2Psca']
for tp in acts_epi_NoSV.obs['timePoint'].cat.categories:
    karthaus_plotting(acts_epi_NoSV,tp,Experiment+'ORA','TypeDiffAmbig','GeneScore')


acts_epi_NoSV.write(f'{Experiment}{epi_ext[:-5]}_Subtypes_ORA.h5ad')





































#epi_gset_dict = {,'Luminal':["Cd24a", "Krt8", "Krt18"]}

# Clean up (removes 1-4 genes per module)
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata_epi.var_names] for key, value in epi_gset_dict.items()})
#print('L3 not detectable')
epi_gset_dict = {key: [item for item in value if item in adata_epi.var_names] for key, value in epi_gset_dict.items()}

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata_epi, gs, ctrl_size=50, n_bins=25, score_name=gset_name)

hist_gene_sigs(adata_epi,epi_gset_dict.keys(),'Epi_Subtypes_ScoreGene',Experiment=Experiment,thresh_gs=None,log_bool=False)

allEpiscores = np.array([adata_epi.obs[gse].tolist() for gse in epi_gset_dict.keys()])
if 'ENZ' in Experiment:
    vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)-1
elif 'SING' in Experiment:
    vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)-1
    vminEpi,vmaxEpi = -0.2,0.6
sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=epi_gset_dict.keys(),vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_Epi_NoSV_Subtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_epi,keys=epi_gset_dict.keys(), groupby=leiRes,save=f'{Experiment}_Epi_ScoreGene_{leiRes}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

MedClusBasal = np.array([np.median(adata_epi.obs['Basal'][adata_epi.obs[leiRes]==clus]) for clus in adata_epi.obs[leiRes].cat.categories])
MedClusLuminal = np.array([np.median(adata_epi.obs['Luminal'][adata_epi.obs[leiRes]==clus]) for clus in adata_epi.obs[leiRes].cat.categories])
#TypeClusters = np.where(MedClusBasal > MedClusLuminal, 'Basal', 'Luminal')
LuminalClusters = np.where(MedClusBasal < MedClusLuminal)[0].astype(str)

adata_epi.obs['EpiLevelLineages'] = 'Basal'
adata_epi.obs.loc[[adata_epi.obs[leiRes].isin(LuminalClusters)][0],'EpiLevelLineages'] = 'Luminal'

sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color='EpiLevelLineages',save=f'{Experiment}_Main_{leiRes}_EpiLevelLineages.png',title=f'{Experiment}_{leiRes}_Epi_Level_Lineages',show=False)

#Subset adata_epi for Luminal clusters...
adata_lum = adata_epi[adata_epi.obs[leiRes].isin(LuminalClusters)].copy()


basal_marker_names = ['Epi_Basal_1','Epi_Basal_SV','Epi_Luminal_SV']
epi_gset_dict = dict()
for gset_name in basal_marker_names:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata_epi.var_names] for key, value in epi_gset_dict.items()})
epi_gset_dict = {key: [item for item in value if item in adata_epi.var_names] for key, value in epi_gset_dict.items()}

#np.intersect1d(sorted(epi_gset_dict['Epi_Basal_1']),sorted(epi_gset_dict['Epi_Basal_SV']))

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata_epi, gs, ctrl_size=len(gs), n_bins=25, score_name=gset_name,use_raw=False)
hist_gene_sigs(adata_epi,basal_marker_names,'Basal_Subtypes_ScoreGene',Experiment=Experiment,thresh_gs=None,log_bool=True)

allEpiscores = np.array([adata_epi.obs[gse].tolist() for gse in basal_marker_names])
vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)
sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=basal_marker_names,vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_Basal_Subtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_epi,keys=basal_marker_names, groupby=leiRes_Prec,save=f'{Experiment}_Basal_ScoreGene_{leiRes_Prec}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

luminal_marker_names = ['Epi_Luminal_1','Epi_Luminal_2Psca']
epi_gset_dict = dict()
for gset_name in luminal_marker_names:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
print('Absent from measured genes:\n',{key: [item for item in value if item not in adata_epi.var_names] for key, value in epi_gset_dict.items()})
epi_gset_dict = {key: [item for item in value if item in adata_epi.var_names] for key, value in epi_gset_dict.items()}

#np.intersect1d(sorted(epi_gset_dict['Epi_Basal_1']),sorted(epi_gset_dict['Epi_Basal_SV']))

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata_epi, gs, ctrl_size=len(gs), n_bins=25, score_name=gset_name,use_raw=False)
hist_gene_sigs(adata_epi,luminal_marker_names,'Luminal_Subtypes_ScoreGene',Experiment=Experiment,thresh_gs=None,log_bool=True)

allEpiscores = np.array([adata_epi.obs[gse].tolist() for gse in luminal_marker_names])
vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)
sc.pl.umap(adata_epi,size=400000/adata_epi.n_obs,alpha=0.8,color=luminal_marker_names,vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_Luminal_Subtypes_ScoreGene.png', cmap='RdBu_r',show=False)

sc.set_figure_params(figsize=(12, 6), dpi_save=300)
sc.pl.violin(adata_epi,keys=luminal_marker_names, groupby=leiRes_Prec,save=f'{Experiment}_Luminal_ScoreGene_{leiRes_Prec}_violin.png',rotation=90,show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)


# Same thing with ORA, the clusters agree....
#adata_epi_ora  = adata_epi.copy()
#adata_epi_ora.obs.drop(columns=['leiden_res1', 'Epi_Basal_SV', 'Epi_Luminal_SV'], inplace=True)
#
#markers_long_list = []
#for gset_name in basal_marker_names:
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
#sc.pl.umap(acts_epi, color=basal_marker_names,save=f'{Experiment}_EpiSV_Subtypes_ORA.png',size=300000/acts_epi.n_obs, cmap='RdBu_r',show=False,vmin=0,vmax=20)
#
#sc.set_figure_params(figsize=(12, 6), dpi_save=300)
#sc.pl.violin(acts_epi,keys=basal_marker_names, groupby=leiRes_Prec,save=f'{Experiment}_SV_ORA_{leiRes_Prec}_violin.png',rotation=90,show=False)
#sc.set_figure_params(figsize=(6, 6), dpi_save=300)
#
#for subtype in basal_marker_names:
#    acts_epi.obs[subtype] = acts_epi.obsm['ora_estimate'][subtype].copy()
#hist_gene_sigs(acts_epi,basal_marker_names,'EpiSV_Subtypes_ORA',Experiment=Experiment,thresh_gs=None,log_bool=True)
#get_highscore_clusters(acts_epi,'Epi_Luminal_SV',leires=leiRes_Prec,medianthresh=median_threshold_SV)
#get_highscore_clusters(acts_epi,'Epi_Basal_SV',leires=leiRes_Prec,medianthresh=median_threshold_SV)

def get_highscore_clusters2(adata_obj,gene_obs,leires='leiden_res2.5',medianthresh=0):
    HighScoreClusters = {}
    for clus in adata_obj.obs[leires].cat.categories:
        MedClus = np.median(adata_obj.obs[gene_obs][adata_obj.obs[leires]==clus])
        print(f'Median Score {gene_obs} {clus}: {MedClus}')
        HighScoreClusters[clus] = MedClus
    return HighScoreClusters

Basal1Clusters = get_highscore_clusters2(adata_epi,'Epi_Basal_1',leires=leiRes_Prec)
BasalSVClusters = get_highscore_clusters2(adata_epi,'Epi_Basal_SV',leires=leiRes_Prec)
LuminalSVClusters = get_highscore_clusters2(adata_epi,'Epi_Luminal_SV',leires=leiRes_Prec)

M = list(Basal1Clusters.values())
Mcl = list(Basal1Clusters.keys())
HighBasal1Clusters = list(np.array(Mcl)[(np.median(M) + 4 * median_abs_deviation(M) < M)])
M = list(BasalSVClusters.values())
Mcl = list(BasalSVClusters.keys())
HighBasalSVClusters = list(np.array(Mcl)[(np.median(M) + 4 * median_abs_deviation(M) < M)])
M = list(LuminalSVClusters.values())
Mcl = list(LuminalSVClusters.keys())
HighLuminalSVClusters = list(np.array(Mcl)[(np.median(M) + 4 * median_abs_deviation(M) < M)])

adata_epi_NoSV = adata_epi[~adata_epi.obs[leiRes_Prec].isin(HighBasal1Clusters+HighBasalSVClusters+HighLuminalSVClusters)].copy()





#How many cells are Epi No SV
#HIST cells in samples, epi vs total 
if 'SING' in Experiment:
    plt.hist(adata.obs['batch'], bins=len(adata.obs['batch'].unique()), alpha=0.5, label='All')#, edgecolor='black')
    plt.hist(adata_epi_NoSV.obs['batch'], bins=len(adata_epi_NoSV.obs['batch'].unique()), alpha=0.5, label='Epithelial')#, edgecolor='black', linewidth=0.1)
elif 'ENZ' in Experiment:
    plt.hist(adata.obs['batch'], bins=len(adata.obs['batch'].unique()), alpha=0.5, label='All', edgecolor='black')
    plt.hist(adata_epi_NoSV.obs['batch'], bins=len(adata_epi_NoSV.obs['batch'].unique()), alpha=0.5, label='Epithelial', edgecolor='black', linewidth=0.3)  
plt.xlabel(None)
plt.ylabel('Number of cells')
plt.title(f'Epithelial cells in {Experiment} Samples')
plt.legend()
plt.xticks(np.arange(0,len(adata_epi_NoSV.obs['batch'].unique()),1)+.5,rotation=45, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.grid(False)
plt.savefig(f'{Experiment}_Main_{leiRes}_ProportionEpiNoSV_Barplot.pdf')
plt.close()

sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["batch"],save=Experiment+'_EPI_NoSV_batch',title=Experiment,show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["timePoint"],save=Experiment+'_EPI_NoSV_timePoint',title=Experiment,show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["tissueProv"],save=Experiment+'_EPI_NoSV_tissueProv',title=Experiment,show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["total_counts", "log1p_total_counts", "n_genes_by_counts", "log1p_n_genes_by_counts"],save=Experiment+'_EPI_NoSV_Counts_QC',show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],save=Experiment+'_EPI_NoSV_PctCounts_QC',show=False)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=["scDblFinder_score","doublet_class"],save=Experiment+'_EPI_NoSV_DbltClass_QC',show=False)

# Writing Epi subset file
adata_epi_NoSV.write(f'{Experiment}{epi_ext}')



### Finer Cell Typing (Basal, L1, L2 esp)

######################
######################
##score_genes METHOD##
######################
######################
##Method 1 with sc.tl.score_genes

adata_epi_NoSV = sc.read_h5ad(f'{Experiment}{epi_ext}')

epi_marker_names_NoSV = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Luminal_3Foxi1']

epi_gset_dict = dict()
for gset_name in epi_marker_names_NoSV:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
epi_gset_dict = {key: [item for item in value if item in adata_epi_NoSV.var_names] for key, value in epi_gset_dict.items()}

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'Gene Scoring {gset_name} ({len(gs)} genes)')
    sc.tl.score_genes(adata_epi_NoSV, gs, ctrl_size=len(gs), n_bins=25, score_name=gset_name)

hist_gene_sigs(adata_epi_NoSV,epi_marker_names_NoSV,'Epi_Subtypes_ScoreGene',Experiment=Experiment,thresh_gs=None,log_bool=False)

allEpiscores = np.array([adata_epi_NoSV.obs[gse].tolist() for gse in epi_marker_names_NoSV])
if 'ENZ' in Experiment:
    vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)-0.4
elif 'SING' in Experiment:
    vminEpi,vmaxEpi = np.min(allEpiscores),np.max(allEpiscores)-1
    vminEpi,vmaxEpi = -0.2,0.6
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=epi_marker_names_NoSV,vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_Epi_NoSV_Subtypes_ScoreGene.png', cmap='RdBu_r',show=False)


#Typing
diff_thresh = 0.05
epi_marker_names_NoSV = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']
# DIFF
typedeflist = []
df_work = adata_epi_NoSV.obs[epi_marker_names_NoSV]
for row in df_work.iterrows():
    markers = ['Basal','L1','L2']

    scores = row[1:][0].tolist()
    maxInd = np.argmax(scores)

    maxType = markers.pop(maxInd)
    maxScore = scores.pop(maxInd)
    AmbigWith = [typ for i,typ in enumerate(markers) if maxScore - diff_thresh < scores[i]]
    
    if AmbigWith:
        typedeflist.append(f"Ambig_{'|'.join(sorted(AmbigWith+[maxType]))}")
    else:
        typedeflist.append(maxType)

adata_epi_NoSV.obs[f'TypeDiffAmbig'] = pd.Categorical(typedeflist,categories=sorted(np.unique(typedeflist), key=len))

#def assign_size(celltype):
#    return 220000/adata_epi_NoSV.n_obs if celltype in ['Basal','L1','L2','L3'] else 120000/adata_epi_NoSV.n_obs
#size_list = adata_epi_NoSV.obs['TypeDiffAmbig'].apply(assign_size).tolist()


adata_epi_NoSV.obs['TypeDiffAmbig'] = quant_singleSeries(adata_epi_NoSV.obs['TypeDiffAmbig'])

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(adata_epi_NoSV,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color=f'TypeDiffAmbig',title=f'{Experiment} ScoreGene Epi Subtypes Min Diff Ambig+',save=f'{Experiment}_Epi_Subtypes_ScoreGene_DiffAmbig.png',show=False)
#,size=size_list

for tp in adata_epi_NoSV.obs['timePoint'].cat.categories:
    karthaus_plotting(adata_epi_NoSV,tp,Experiment,'TypeDiffAmbig','GeneScore')

adata_epi_NoSV.write(f'{Experiment}{epi_ext[:-5]}_Subtypes_ScoreGene.h5ad')


#def assign_ambig(celltype):
#    return 'Ambig' if 'Ambig_' in celltype else celltype
#
#less_ambig = adata_epi_NoSV.obs['TypeDiffAmbig'].apply(assign_ambig)
#adata_epi_NoSV.obs['TypeDiff'] = pd.Categorical(less_ambig, categories=sorted(np.unique(less_ambig), key=len))
#sc.pl.umap(adata_epi_NoSV,size=220000/adata_epi_NoSV.n_obs,color=f'TypeDiff',title=f'{Experiment} ScoreGene Epi Subtypes Min Diff',save=f'{Experiment}_Epi_Subtypes_ScoreGene_Diff.png',show=False)
#
#
## ITER
##print(f'\nSeparating Untyped')
#step_thresh,max_thresh,min_thresh=0.02,2,0
#base_thresh = Subtypes_Thresh_ScatterL1L2B(adata_epi_NoSV,'Epi_NoSV_ScoreGene',Experiment=Experiment)
#ada_typ = adata_epi_NoSV.copy()
#aaTypes = {'L1':[],'L2':[],'Basal':[],'Ambig':[]}
#new_thresh = base_thresh
#Untyped={0,1}
#while (new_thresh > min_thresh) & (len(Untyped)>0):
#    #print(f'Threshold={new_thresh}')
#    
#    setL1 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_1'] >= new_thresh])
#    setL2 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_2Psca'] >= new_thresh])
#    setb = set(ada_typ.obs.index[ada_typ.obs['Epi_Basal_1'] >= new_thresh])
#
#    Ambig = set(setb & setL1 | setb & setL2 | setL1 & setL2)
#    Typed = set(setL1 | setL2 | setb)
#
#    aaTypes['L1'].extend(list(setL1 - Ambig))
#    aaTypes['L2'].extend(list(setL2 - Ambig))
#    aaTypes['Basal'].extend(list(setb - Ambig))
#    aaTypes['Ambig'].extend(list(Ambig))
#
#    Untyped = set(ada_typ.obs.index) - Typed
#    aaTypes['Untyped'] = list(Untyped)
#
#    tot_cel_typ = sum([len(aaTypes[k]) for k in aaTypes.keys()])
#
#    assert tot_cel_typ==adata_epi_NoSV.n_obs,'Something is wrong with the type separation'
#
#    ada_typ = ada_typ[ada_typ.obs.index.isin(Untyped)].copy()
#
#    new_thresh -= step_thresh
#    #print(len(Untyped))
#Untyped_LowScore = set(aaTypes['Untyped'])
#
##print(f'\nResult:')
##for k in aaTypes.keys():
##    print(f'{k}: {len(aaTypes[k])}')
##
##print(f'\nSeparating Ambigs')
#ada_typ = adata_epi_NoSV[adata_epi_NoSV.obs.index.isin(aaTypes['Ambig'])].copy()
#new_thresh = base_thresh
#while (len(aaTypes['Ambig']) > 0) & (new_thresh <= max_thresh):
#    ambigbefore = len(aaTypes['Ambig'])
#    #print(f'Threshold={new_thresh}')
#    
#    setL1 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_1'] >= new_thresh])
#    setL2 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_2Psca'] >= new_thresh])
#    setb = set(ada_typ.obs.index[ada_typ.obs['Epi_Basal_1'] >= new_thresh])
#
#    Typed = set(setL1 | setL2 | setb)
#    Ambig = set(setb & setL1 | setb & setL2 | setL1 & setL2)
#    aaTypes['Ambig'] = list(Ambig)
#
#    aaTypes['L1'].extend(list(setL1 - Ambig))
#    aaTypes['L2'].extend(list(setL2 - Ambig))
#    aaTypes['Basal'].extend(list(setb - Ambig))
#
#    Untyped = set(ada_typ.obs.index) - Typed
#    aaTypes['Untyped'].extend(list(Untyped))
#
#    tot_cel_typ = sum([len(aaTypes[k]) for k in aaTypes.keys()])
#
#    assert tot_cel_typ==adata_epi_NoSV.n_obs,'Something is wrong with the type separation'
#
#    ada_typ = ada_typ[ada_typ.obs.index.isin(Ambig)].copy()
#
#    new_thresh += step_thresh
#
#    diffAmbig = ambigbefore - len(aaTypes['Ambig'])
#
#    #print(diffAmbig)
#
#if len(aaTypes['Ambig'])==0:
#    aaTypes.pop('Ambig')
#
##print(f'\nResult:')
##for k in aaTypes.keys():
##    print(f'{k}: {len(aaTypes[k])}')
#
#aaTypes['Untyped_Ambig'] = list(set(aaTypes.pop('Untyped')) - Untyped_LowScore)
#aaTypes['Untyped_LowScore'] = list(Untyped_LowScore)
#
#adata_epi_NoSV.obs[f'TypeIterAmbig'] = pd.Categorical([list(aaTypes.keys())[0]]*adata_epi_NoSV.n_obs, categories=sorted(np.unique(list(aaTypes.keys())), key=len))
#for gset_name,TypeDef_cells_gset_name in aaTypes.items():
#    adata_epi_NoSV.obs.loc[adata_epi_NoSV.obs.index.isin(TypeDef_cells_gset_name), 'TypeIterAmbig'] = gset_name
#
#def assign_size(celltype):
#    return 220000/adata_epi_NoSV.n_obs if celltype in ['Basal','L1','L2','L3'] else 120000/adata_epi_NoSV.n_obs
#size_list = adata_epi_NoSV.obs['TypeIterAmbig'].apply(assign_size).tolist()
#
#sc.pl.umap(adata_epi_NoSV,size=size_list,color=f'TypeIterAmbig',title=f'{Experiment} ScoreGene Epi Subtypes Iter Thresh Ambig+',save=f'{Experiment}_Epi_Subtypes_ScoreGene_IterAmbig.png',show=False)
#
#aaTypes['Ambig'] = aaTypes.pop('Untyped_Ambig') + aaTypes.pop('Untyped_LowScore')
#
#adata_epi_NoSV.obs[f'TypeIter'] = pd.Categorical([list(aaTypes.keys())[0]]*adata_epi_NoSV.n_obs, categories=sorted(np.unique(list(aaTypes.keys())), key=len))
#for gset_name,TypeDef_cells_gset_name in aaTypes.items():
#    adata_epi_NoSV.obs.loc[adata_epi_NoSV.obs.index.isin(TypeDef_cells_gset_name), 'TypeIter'] = gset_name
#sc.pl.umap(adata_epi_NoSV,size=220000/adata_epi_NoSV.n_obs,color=f'TypeIter',title=f'{Experiment} ScoreGene Epi Subtypes Iter Thresh',save=f'{Experiment}_Epi_Subtypes_ScoreGene_Iter.png',show=False)
#
#
#adata_epi_NoSV.write(f'{Experiment}{epi_ext[:-5]}_Subtypes_ScoreGene.h5ad')



######################
######################
######ORA METHOD######
######################
######################
#Method 2 with dc.run_ora

#EPI
adata_epi_NoSV = sc.read_h5ad(f'{Experiment}{epi_ext}')

epi_marker_names_NoSV = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Luminal_3Foxi1']

epi_gset_dict = dict()
for gset_name in epi_marker_names_NoSV:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
epi_gset_dict = {key: [item for item in value if item in adata_epi_NoSV.var_names] for key, value in epi_gset_dict.items()}

markers_long_list = []
for gset_name in epi_marker_names_NoSV:
    for gen in list(gset_df[gset_name].dropna()):
        markers_long_list.append([gen,gset_name])

markers_long = pd.DataFrame(markers_long_list,columns = ['genesymbol', 'cell_type'])

dc.run_ora(
    mat=adata_epi_NoSV,
    net=markers_long,
    source='cell_type',
    target='genesymbol',
    min_n=10,
    verbose=True,
    use_raw=False
)

acts_epi_NoSV = dc.get_acts(adata_epi_NoSV, obsm_key='ora_estimate')
# We need to remove inf and set them to the maximum value observed for pvals=0
acts_epi_v = acts_epi_NoSV.X.ravel()
max_e = np.nanmax(acts_epi_v[np.isfinite(acts_epi_v)])
acts_epi_NoSV.X[~np.isfinite(acts_epi_NoSV.X)] = max_e

if 'ENZ' in Experiment:
    vminEpi,vmaxEpi = 0,30
elif 'SING' in Experiment:
    vminEpi,vmaxEpi = 0,15

sc.pl.umap(acts_epi_NoSV, color=epi_marker_names_NoSV,save=f'{Experiment}_Epi_NoSV_Subtypes_ORA.png',size=300000/acts_epi_NoSV.n_obs, cmap='RdBu_r',show=False,vmin=vminEpi,vmax=vmaxEpi)

for subtype in epi_marker_names_NoSV:
    acts_epi_NoSV.obs[subtype] = acts_epi_NoSV.obsm['ora_estimate'][subtype].copy()
hist_gene_sigs(acts_epi_NoSV,epi_marker_names_NoSV,'Epi_Subtypes_ORA',Experiment=Experiment,thresh_gs=None,log_bool=True)

#Typing
diff_thresh = 3
step_thresh,max_thresh,min_thresh=1,50,4
epi_marker_names_NoSV = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']

# DIFF
typedeflist = []
df_work = acts_epi_NoSV.obs[epi_marker_names_NoSV]
for row in df_work.iterrows():
    markers = ['Basal','L1','L2']

    scores = row[1:][0].tolist()
    maxInd = np.argmax(scores)

    maxType = markers.pop(maxInd)
    maxScore = scores.pop(maxInd)
    AmbigWith = [typ for i,typ in enumerate(markers) if maxScore - diff_thresh < scores[i]]
    
    if AmbigWith:
        typedeflist.append(f"Ambig_{'|'.join(sorted(AmbigWith+[maxType]))}")
    else:
        typedeflist.append(maxType)

acts_epi_NoSV.obs[f'TypeDiffAmbig'] = pd.Categorical(typedeflist,categories=sorted(np.unique(typedeflist), key=len))

acts_epi_NoSV.obs['TypeDiffAmbig'] = quant_singleSeries(acts_epi_NoSV.obs['TypeDiffAmbig'])

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(acts_epi_NoSV,size=400000/acts_epi_NoSV.n_obs,alpha=0.8,color=f'TypeDiffAmbig',title=f'{Experiment} ScoreGene Epi Subtypes Min Diff Ambig+',save=f'{Experiment}_Epi_Subtypes_ScoreGene_DiffAmbig.png',show=False)
#,size=size_list

for tp in adata_epi_NoSV.obs['timePoint'].cat.categories:
    karthaus_plotting(adata_epi_NoSV,tp,Experiment,'TypeDiffAmbig','GeneScore')


acts_epi_NoSV.write(f'{Experiment}{epi_ext[:-5]}_Subtypes_ORA.h5ad')










def assign_size(celltype):
    return 220000/acts_epi_NoSV.n_obs if celltype in ['Basal','L1','L2','L3'] else 120000/acts_epi_NoSV.n_obs
size_list = acts_epi_NoSV.obs['TypeDiffAmbig'].apply(assign_size).tolist()

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
sc.pl.umap(acts_epi_NoSV,color=f'TypeDiffAmbig',title=f'{Experiment} ORA Epi Subtypes Min Diff Ambig+',save=f'{Experiment}_Epi_Subtypes_ORA_DiffAmbig.png',show=False)
#,size=size_list

def assign_ambig(celltype):
    return 'Ambig' if 'Ambig_' in celltype else celltype

less_ambig = acts_epi_NoSV.obs['TypeDiffAmbig'].apply(assign_ambig)
acts_epi_NoSV.obs['TypeDiff'] = pd.Categorical(less_ambig, categories=sorted(np.unique(less_ambig), key=len))
sc.pl.umap(acts_epi_NoSV,size=220000/acts_epi_NoSV.n_obs,color=f'TypeDiff',title=f'{Experiment} ORA Epi Subtypes Min Diff',save=f'{Experiment}_Epi_Subtypes_ORA_Diff.png',show=False)


# ITER
# Cant do it in SING because the starting score is too low
if 'ENZ' in Experiment:
    print(f'\nSeparating Untyped')
    base_thresh = Subtypes_Thresh_ScatterL1L2B(acts_epi_NoSV,'Epi_NoSV_ORA',Experiment=Experiment,thresholds = [i for i in range(0, 50)])
    ada_typ = acts_epi_NoSV.copy()
    aaTypes = {'L1':[],'L2':[],'Basal':[],'Ambig':[]}
    new_thresh = base_thresh
    Untyped={0,1}
    while (new_thresh > min_thresh) & (len(Untyped)>0):
        #print(f'Threshold={new_thresh}')
        
        setL1 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_1'] >= new_thresh])
        setL2 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_2Psca'] >= new_thresh])
        setb = set(ada_typ.obs.index[ada_typ.obs['Epi_Basal_1'] >= new_thresh])

        Ambig = set(setb & setL1 | setb & setL2 | setL1 & setL2)
        Typed = set(setL1 | setL2 | setb)

        aaTypes['L1'].extend(list(setL1 - Ambig))
        aaTypes['L2'].extend(list(setL2 - Ambig))
        aaTypes['Basal'].extend(list(setb - Ambig))
        aaTypes['Ambig'].extend(list(Ambig))

        Untyped = set(ada_typ.obs.index) - Typed
        aaTypes['Untyped'] = list(Untyped)

        tot_cel_typ = sum([len(aaTypes[k]) for k in aaTypes.keys()])

        assert tot_cel_typ==acts_epi_NoSV.n_obs,'Something is wrong with the type separation'

        ada_typ = ada_typ[ada_typ.obs.index.isin(Untyped)].copy()

        new_thresh -= step_thresh
        #print(len(Untyped))
    Untyped_LowScore = set(aaTypes['Untyped'])

    #print(f'\nResult:')
    #for k in aaTypes.keys():
    #    print(f'{k}: {len(aaTypes[k])}')
    #
    #print(f'\nSeparating Ambigs')
    ada_typ = acts_epi_NoSV[acts_epi_NoSV.obs.index.isin(aaTypes['Ambig'])].copy()
    new_thresh = base_thresh
    while (len(aaTypes['Ambig']) > 0) & (new_thresh <= max_thresh):
        ambigbefore = len(aaTypes['Ambig'])
        #print(f'Threshold={new_thresh}')
        
        setL1 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_1'] >= new_thresh])
        setL2 = set(ada_typ.obs.index[ada_typ.obs['Epi_Luminal_2Psca'] >= new_thresh])
        setb = set(ada_typ.obs.index[ada_typ.obs['Epi_Basal_1'] >= new_thresh])

        Typed = set(setL1 | setL2 | setb)
        Ambig = set(setb & setL1 | setb & setL2 | setL1 & setL2)
        aaTypes['Ambig'] = list(Ambig)

        aaTypes['L1'].extend(list(setL1 - Ambig))
        aaTypes['L2'].extend(list(setL2 - Ambig))
        aaTypes['Basal'].extend(list(setb - Ambig))

        Untyped = set(ada_typ.obs.index) - Typed
        aaTypes['Untyped'].extend(list(Untyped))

        tot_cel_typ = sum([len(aaTypes[k]) for k in aaTypes.keys()])

        assert tot_cel_typ==acts_epi_NoSV.n_obs,'Something is wrong with the type separation'

        ada_typ = ada_typ[ada_typ.obs.index.isin(Ambig)].copy()

        new_thresh += step_thresh

        diffAmbig = ambigbefore - len(aaTypes['Ambig'])

        #print(diffAmbig)

    if len(aaTypes['Ambig'])==0:
        aaTypes.pop('Ambig')

    #print(f'\nResult:')
    #for k in aaTypes.keys():
    #    print(f'{k}: {len(aaTypes[k])}')

    aaTypes['Untyped_Ambig'] = list(set(aaTypes.pop('Untyped')) - Untyped_LowScore)
    aaTypes['Untyped_LowScore'] = list(Untyped_LowScore)

    acts_epi_NoSV.obs[f'TypeIterAmbig'] = pd.Categorical([list(aaTypes.keys())[0]]*acts_epi_NoSV.n_obs, categories=sorted(np.unique(list(aaTypes.keys())), key=len))
    for gset_name,TypeDef_cells_gset_name in aaTypes.items():
        acts_epi_NoSV.obs.loc[acts_epi_NoSV.obs.index.isin(TypeDef_cells_gset_name), 'TypeIterAmbig'] = gset_name

    def assign_size(celltype):
        return 220000/acts_epi_NoSV.n_obs if celltype in ['Basal','L1','L2','L3'] else 120000/acts_epi_NoSV.n_obs
    size_list = acts_epi_NoSV.obs['TypeIterAmbig'].apply(assign_size).tolist()

    sc.pl.umap(acts_epi_NoSV,size=size_list,color=f'TypeIterAmbig',title=f'{Experiment} ORA Epi Subtypes Iter Thresh Ambig+',save=f'{Experiment}_Epi_Subtypes_ORA_IterAmbig.png',show=False)

    aaTypes['Ambig'] = aaTypes.pop('Untyped_Ambig') + aaTypes.pop('Untyped_LowScore')

    acts_epi_NoSV.obs[f'TypeIter'] = pd.Categorical([list(aaTypes.keys())[0]]*acts_epi_NoSV.n_obs, categories=sorted(np.unique(list(aaTypes.keys())), key=len))
    for gset_name,TypeDef_cells_gset_name in aaTypes.items():
        acts_epi_NoSV.obs.loc[acts_epi_NoSV.obs.index.isin(TypeDef_cells_gset_name), 'TypeIter'] = gset_name
    sc.pl.umap(acts_epi_NoSV,size=220000/acts_epi_NoSV.n_obs,color=f'TypeIter',title=f'{Experiment} ORA Epi Subtypes Iter Thresh',save=f'{Experiment}_Epi_Subtypes_ORA_Iter.png',show=False)


acts_epi_NoSV.write(f'{Experiment}{epi_ext[:-5]}_Subtypes_ORA.h5ad')


######################
######################
####METHOD OVERLAP####
######################
######################


### CellTyping Overlap analysis

ora = sc.read_h5ad(f'{Experiment}{epi_ext[:-5]}_Subtypes_ORA.h5ad')
genescore = sc.read_h5ad(f'{Experiment}{epi_ext[:-5]}_Subtypes_ScoreGene.h5ad')

for gset_name in ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']:
    print(f'Plotting {gset_name} ORA_GeneScore Scatter')

    oradf = pd.DataFrame({'obs_names': ora.obs_names, gset_name:ora.obs[gset_name], 'celltype':ora.obs['TypeDiffAmbig']})
    genescoredf = pd.DataFrame({'obs_names': genescore.obs_names, gset_name:genescore.obs[gset_name], 'celltype':genescore.obs['TypeDiffAmbig']})
    merged_df = oradf.merge(genescoredf, on='obs_names', how='inner')


    if gset_name == 'Epi_Basal_1':
        main = 'Basal'
    elif gset_name == 'Epi_Luminal_1':
        main = 'L1'
    elif gset_name == 'Epi_Luminal_2Psca':
        main = 'L2'

    co = []
    for cx,cy in zip(merged_df['celltype_x'],merged_df['celltype_y']):
        if (main in cx) | (main in cy):
            if cx==cy:
                co.append(cx)
            else:
                co.append(f'{main} discordance')
        else:
            co.append(f'Other type')

    # Extract unique categories from co
    unique_categories = np.unique(co)

    # Generate a colormap with a color for each unique category
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_categories)))

    # Create a dictionary to map each category to its color
    category_color_map = dict(zip(unique_categories, colors))

    # Plot points for each category separately with its own color
    for category in unique_categories:
        category_indices = [i for i, val in enumerate(co) if val == category]
        plt.scatter(merged_df.iloc[category_indices][f'{gset_name}_x'], 
                    merged_df.iloc[category_indices][f'{gset_name}_y'], 
                    s=6, color=category_color_map[category],alpha=0.4, label=category)

    # Label axes and set a title
    plt.xlabel(f'ORA Score')
    plt.ylabel(f'Gene Score')
    plt.title(f'{gset_name}')
    plt.suptitle(f'{Experiment}')
    plt.legend(fontsize='small')

    
    # Display the Pearson correlation coefficient in a box in the top right corner
    plt.text(0.85, 0.9, f'r = {merged_df[f'{gset_name}_x'].corr(merged_df[f'{gset_name}_y']):.2f}', transform=plt.gca().transAxes,
        bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='round,pad=0.3'), fontsize=10)

    # Set x and y-axis limits
    #plt.xlim(rnaxlim[0], rnaxlim[1])  # Adjust the limits based on data
    #plt.ylim(atacylim[0], atacylim[1])  # Adjust the limits based on data

    # Display the plot
    plt.tight_layout()

    # Save the scatter plot
    plt.savefig(f'{Experiment}_ORA_GeneScore_{gset_name}_Scatter.pdf')
    plt.close()


#########
#Matrix comparison Scoring Method
#########

for TypeMethod in ['TypeDiff']: #TypeIter
    print(f'Plotting {Experiment} {TypeMethod}')

    oradf = pd.DataFrame({'Experiment': Experiment, 'obs_names': ora.obs_names, TypeMethod: ora.obs[TypeMethod], 'Scoring':f'ORA_{TypeMethod}'})
    #oradf[TypeMethod] = oradf[TypeMethod].str.split(r' \(').str.get(0)
    genescoredf = pd.DataFrame({'Experiment': Experiment, 'obs_names': genescore.obs_names, TypeMethod: genescore.obs[TypeMethod], 'Scoring':f'ScoreGene_{TypeMethod}'})
    #genescoredf[TypeMethod] = genescoredf[TypeMethod].str.split(r' \(').str.get(0)

    merged_df = oradf.merge(genescoredf, on='obs_names', how='inner')

    groups1 = np.array(merged_df[f'{TypeMethod}_x'])
    cl1 = np.unique(groups1)
    groups2 = np.array(merged_df[f'{TypeMethod}_y'])
    cl2 = np.unique(groups2)

    # create a heatmap of the comparison
    # Using Jaccard Index (inters(A,B)/union(A,B))
    heatmap_data = np.zeros((len(cl1), len(cl2)))

    for i, group1 in enumerate(cl1):
        size1 = len(np.where(groups1==group1)[0])

        for j, group2 in enumerate(cl2):
            size2 = len(np.where(groups2==group2)[0])

            intersection = len(np.intersect1d(np.where(groups1==group1), np.where(groups2==group2)))
            heatmap_data[i, j] = intersection / (size1 + size2 - intersection)

    #Sort x axis
    max_x = np.max(heatmap_data,axis=0)
    sorted_indices_x = np.argsort(-max_x)
    #Sort y axis
    max_y = np.max(heatmap_data,axis=1)
    sorted_indices_y = np.argsort(-max_y)

    #Sort heatmap
    heatmap_data = heatmap_data[sorted_indices_y][:, sorted_indices_x]
     
    #plt.figure(figsize=(10, 6),dpi=300)
    fig, ax = plt.subplots(figsize=(10, 6),dpi=300)
    im = ax.imshow(heatmap_data, cmap='coolwarm', vmin=0, vmax=1)

    # set up the plot
    ax.set_xticks(np.arange(len(cl2)))
    ax.set_yticks(np.arange(len(cl1)))
    #Put ticks in right order 
    ax.set_xticklabels(cl2[sorted_indices_x])
    ax.set_yticklabels(cl1[sorted_indices_y])

    #rotate xticks
    ax.tick_params(axis='x', rotation=90)

    # Add the number of observations in each group on x and y axes
    for i, (value, label) in enumerate(zip(cl2[sorted_indices_x], ax.get_xticklabels())):
        ax.text(i, -0.56, f'{len(np.where(groups2 == value)[0])}', ha='left', va='center',rotation=45,rotation_mode='anchor')

    for i, (value, label) in enumerate(zip(cl1[sorted_indices_y], ax.get_yticklabels())):
        ax.text(-0.45, i, f'{len(np.where(groups1 == value)[0])}', ha='left', va='center')
        label.set_y(i - 0.3)  # Adjust the position of the current label

    ax.set_xlabel('GeneScore Epithelial Subtypes')
    ax.set_ylabel('ORA Epithelial Subtypes')
    ax.set_title(f'{Experiment} {TypeMethod}',pad=40)
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.grid(False)
    plt.tight_layout()

    # save the plot
    plt.savefig(f'{Experiment}_ORA_GeneScore_Jaccard_{TypeMethod}.pdf')
    plt.close()


#########
#Matrix comparison Typing Method
#########

for ScoringMethod in ['ScoreGene','ORA']:
    if ('SING' in Experiment) & (ScoringMethod=='ORA'):
        continue
    
    path_adata = f'{Experiment}{epi_ext[:-5]}_Subtypes_{ScoringMethod}.h5ad'
    adat = sc.read_h5ad(f'./{path_adata}')
    print(f'Plotting {Experiment} {ScoringMethod}')
    TypeMethod1 = 'TypeDiff'    
    TypeMethod2 = 'TypeIter'


    df1 = pd.DataFrame({'Experiment': Experiment, 'obs_names': adat.obs_names, TypeMethod1: adat.obs[TypeMethod1], 'Scoring':f'{TypeMethod1}'})
    df2 = pd.DataFrame({'Experiment': Experiment, 'obs_names': adat.obs_names, TypeMethod2: adat.obs[TypeMethod2], 'Scoring':f'{TypeMethod2}'})
    merged_df = df1.merge(df2, on='obs_names', how='inner')

    groups1 = np.array(merged_df[f'{TypeMethod1}'])
    cl1 = np.unique(groups1)
    groups2 = np.array(merged_df[f'{TypeMethod2}'])
    cl2 = np.unique(groups2)

    # create a heatmap of the comparison
    # Using Jaccard Index (inters(A,B)/union(A,B))
    heatmap_data = np.zeros((len(cl1), len(cl2)))

    for i, group1 in enumerate(cl1):
        size1 = len(np.where(groups1==group1)[0])

        for j, group2 in enumerate(cl2):
            size2 = len(np.where(groups2==group2)[0])

            intersection = len(np.intersect1d(np.where(groups1==group1), np.where(groups2==group2)))
            heatmap_data[i, j] = intersection / (size1 + size2 - intersection)

    #Sort x axis
    max_x = np.max(heatmap_data,axis=0)
    sorted_indices_x = np.argsort(-max_x)
    #Sort y axis
    max_y = np.max(heatmap_data,axis=1)
    sorted_indices_y = np.argsort(-max_y)

    #Sort heatmap
    heatmap_data = heatmap_data[sorted_indices_y][:, sorted_indices_x]
     
    #plt.figure(figsize=(10, 6),dpi=300)
    fig, ax = plt.subplots(figsize=(10, 6),dpi=300)
    im = ax.imshow(heatmap_data, cmap='coolwarm', vmin=0, vmax=1)

    # set up the plot
    ax.set_xticks(np.arange(len(cl2)))
    ax.set_yticks(np.arange(len(cl1)))
    #Put ticks in right order 
    ax.set_xticklabels(cl2[sorted_indices_x])
    ax.set_yticklabels(cl1[sorted_indices_y])

    #rotate xticks
    ax.tick_params(axis='x', rotation=90)

    # Add the number of observations in each group on x and y axes
    for i, (value, label) in enumerate(zip(cl2[sorted_indices_x], ax.get_xticklabels())):
        ax.text(i, -0.56, f'{len(np.where(groups2 == value)[0])}', ha='left', va='center',rotation=45,rotation_mode='anchor')

    for i, (value, label) in enumerate(zip(cl1[sorted_indices_y], ax.get_yticklabels())):
        ax.text(-0.45, i, f'{len(np.where(groups1 == value)[0])}', ha='left', va='center')
        label.set_y(i - 0.3)  # Adjust the position of the current label

    ax.set_xlabel(f'{TypeMethod2} Epithelial Subtypes')
    ax.set_ylabel(f'{TypeMethod1} Epithelial Subtypes')
    ax.set_title(f'{Experiment} {ScoringMethod} Diff vs Iter',pad=40)
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.grid(False)
    plt.tight_layout()

    # save the plot
    plt.savefig(f'{Experiment}_{ScoringMethod}_Jaccard_DiffvsIter.pdf')
    plt.close()





































#OLD SCORE GENE

if 'ENZ' in Experiment:
    threshold_genescoreSub = 0.40
elif 'SING' in Experiment:
    threshold_genescoreSub = 0.55

venn3_subtypes(adata_epi_NoSV,epi_marker_names_NoSV,threshold_genescoreSub,Experiment,'GeneScore')

setL1 = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Luminal_1'] > threshold_genescoreSub])
setL2 = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Luminal_2Psca'] > threshold_genescoreSub])
setb = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Basal_1'] > threshold_genescoreSub])

All_Ambig = set(setb & setL1 | setb & setL2 | setL1 & setL2)
L1L2B_Ambig = set(setb & setL1 & setL2)
L1L2_Ambig = set(setL1 & setL2) - L1L2B_Ambig
L1B_Ambig = set(setL1 & setb) - L1L2B_Ambig
L2B_Ambig = set(setb & setL2) - L1L2B_Ambig
setL1_NoAmbig = setL1 - All_Ambig
setL2_NoAmbig = setL2 - All_Ambig
setb_NoAmbig = setb - All_Ambig

Epi_cats = epi_marker_names_NoSV+['L1/L2_Ambig','L1/B_Ambig','L2/B_Ambig','L1/L2/B_Ambig','Untyped']
Epi_cats_sets = [setb_NoAmbig,setL1_NoAmbig,setL2_NoAmbig,L1L2_Ambig,L1B_Ambig,L2B_Ambig,L1L2B_Ambig]
for i,TypeDef_cells_gset_name in enumerate(Epi_cats_sets):
    Epi_cats[i] = Epi_cats[i] + f' ({len(TypeDef_cells_gset_name)} cells)'
Epi_cats[-1] = Epi_cats[-1] + f' ({adata_epi_NoSV.n_obs - len(setL1 | setL2 | setb)} cells)'

adata_epi_NoSV.obs[f'TypeDef'] = pd.Categorical([Epi_cats[-1]]*adata_epi_NoSV.n_obs, categories=Epi_cats)

for gset_name,TypeDef_cells_gset_name in zip(Epi_cats[:-1],Epi_cats_sets):
    adata_epi_NoSV.obs.loc[adata_epi_NoSV.obs.index.isin(TypeDef_cells_gset_name), 'TypeDef'] = gset_name

#TypeDef_cells_gset_name = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs[gset_name] > threshold_genescoreSub]) - setAmbig
size_list = []
for cel in adata_epi_NoSV.obs.index:
    if [item for item in epi_marker_names_NoSV if adata_epi_NoSV.obs['TypeDef'][cel][:10] in item]:
        size_list.append(360000/adata_epi_NoSV.n_obs)
    else:
        size_list.append(220000/adata_epi_NoSV.n_obs)

sc.pl.umap(adata_epi_NoSV,size=size_list,color=f'TypeDef',alpha=0.7,title=f'{Experiment}_Epi_Subtypes',save=f'{Experiment}_Epi_Subtypes_ScoreGene_TypeDef.png',show=False)
sc.pl.umap(adata_epi_NoSV,color=epi_marker_names_NoSV,vmin=vminEpi,vmax=vmaxEpi,save=f'{Experiment}_Epi_Subtypes_NoSV_ScoreGene.png',size=300000/adata_epi_NoSV.n_obs, cmap='RdBu_r',show=False)

#sc.set_figure_params(figsize=(8, 6), dpi_save=300)
#for subtype in epi_marker_names_NoSV:
#    mean_values = adata_epi_NoSV.obs.groupby(leiRes,observed=False)[subtype].median()
#    sorted_clusters = mean_values.sort_values(ascending=False).index.tolist()
#    sc.pl.violin(adata_epi_NoSV, keys=subtype, groupby=leiRes,order=sorted_clusters,save=f'{Experiment}_{subtype}_{leiRes}_violin.png',show=False)
#sc.set_figure_params(figsize=(6, 6), dpi_save=300)

adata_epi_NoSV.write(f'{Experiment}{epi_ext[:-5]}_Subtypes_ScoreGene.h5ad')







# OLD ORA

epi_marker_names_NoSV = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']

if 'ENZ' in Experiment:
    threshold_genescoreSub = threshold_SV = 12
elif 'SING' in Experiment:
    threshold_genescoreSub = 12

venn3_subtypes(acts_epi_NoSV,epi_marker_names_NoSV,threshold_genescoreSub,Experiment,'ORA')

adata_epi_NoSV = acts_epi_NoSV

setL1 = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Luminal_1'] > threshold_genescoreSub])
setL2 = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Luminal_2Psca'] > threshold_genescoreSub])
setb = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Basal_1'] > threshold_genescoreSub])

All_Ambig = set(setb & setL1 | setb & setL2 | setL1 & setL2)
L1L2B_Ambig = set(setb & setL1 & setL2)
L1L2_Ambig = set(setL1 & setL2) -  L1L2B_Ambig
L1B_Ambig = set(setL1 & setb) -  L1L2B_Ambig
L2B_Ambig = set(setb & setL2) -  L1L2B_Ambig
setL1_NoAmbig = setL1 - All_Ambig
setL2_NoAmbig = setL2 - All_Ambig
setb_NoAmbig = setb - All_Ambig

Epi_cats = epi_marker_names_NoSV+['L1/L2_Ambig','L1/B_Ambig','L2/B_Ambig','L1/L2/B_Ambig','Untyped']
Epi_cats_sets = [setb_NoAmbig,setL1_NoAmbig,setL2_NoAmbig,L1L2_Ambig,L1B_Ambig,L2B_Ambig,L1L2B_Ambig]
for i,TypeDef_cells_gset_name in enumerate(Epi_cats_sets):
    Epi_cats[i] = Epi_cats[i] + f' ({len(TypeDef_cells_gset_name)} cells)'
Epi_cats[-1] = Epi_cats[-1] + f' ({adata_epi_NoSV.n_obs - len(setL1 | setL2 | setb)} cells)'

adata_epi_NoSV.obs[f'TypeDef'] = pd.Categorical([Epi_cats[-1]]*adata_epi_NoSV.n_obs, categories=Epi_cats)
for gset_name,TypeDef_cells_gset_name in zip(Epi_cats[:-1],Epi_cats_sets):
    adata_epi_NoSV.obs.loc[adata_epi_NoSV.obs.index.isin(TypeDef_cells_gset_name), 'TypeDef'] = gset_name

#TypeDef_cells_gset_name = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs[gset_name] > threshold_genescoreSub]) - setAmbig
size_list = []
for cel in adata_epi_NoSV.obs.index:
    if [item for item in epi_marker_names_NoSV if adata_epi_NoSV.obs['TypeDef'][cel][:10] in item]:
        size_list.append(360000/adata_epi_NoSV.n_obs)
    else:
        size_list.append(220000/adata_epi_NoSV.n_obs)

sc.set_figure_params(figsize=(6, 6), dpi_save=300)



#SV removal at the single cell level
SV_cells = set(adata_epi.obs.index[adata_epi.obs['Epi_Basal_SV'] >= threshold_SV]) | set(adata_epi.obs.index[adata_epi.obs['Epi_Luminal_SV'] >= threshold_SV])
adata_epi_NoSV = adata_epi[~adata_epi.obs.index.isin(SV_cells)].copy()

threshold_genescoreSub = Subtypes_Thresh_ScatterL1L2B(adata_epi_NoSV,'Epi',Experiment=Experiment)

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

    plt.scatter(x=thresholds,y=np.array(bSV), alpha=0.8, color='#a6cee3', edgecolor='black', linewidth=0.2, label="Basal SV only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(lSV), alpha=0.8, color='#b2df8a', edgecolor='black', linewidth=0.2, label="Luminal SV only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(ambig), alpha=0.8, color='#bdbdbd', edgecolor='black', linewidth=0.2, label="Ambig",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(NoTypes), alpha=0.8, color='#fb9a99', edgecolor='black', linewidth=0.2, label="Untyped",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.yscale('log')
    plt.ylabel(f'% of Epithelial cells')
    plt.xlabel('Score threshold for subtype assignment')
    plt.title('Uniquely assigned Epithelial subtypes')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{Experiment}_{name}_NCells_SV_Scatter.pdf')
    plt.close()

    return thresholds[np.argmin(np.array(ambig)+np.array(NoTypes))]

threshold_SV = Subtypes_Thresh_ScatterSV(adata_epi,'SV',Experiment=Experiment)







# Cluster based type defing with ORA
#leiRes='leiden_res3'
#sc.set_figure_params(figsize=(8, 6), dpi_save=300)
#for subtype in epi_marker_names_NoSV:
#    acts_epi.obs[subtype] = acts_epi.obsm['ora_estimate'][subtype].copy()
#    mean_values = acts_epi.obs.groupby(leiRes,observed=False)[subtype].median()
#    sorted_clusters = mean_values.sort_values(ascending=False).index.tolist()
#    acts_epi.obs.drop(columns=subtype, inplace=True)
#    sc.pl.violin(acts_epi, keys=subtype, groupby=leiRes,order=sorted_clusters,save=f'{Experiment}_{subtype}_{leiRes}_EpiORA_violin.png',show=False)
#sc.set_figure_params(figsize=(6, 6), dpi_save=300)

#df = dc.rank_sources_groups(acts_epi, groupby=leiRes, reference='rest', method='t-test_overestim_var')

#n_ctypes = 3
#ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
#sc.set_figure_params(figsize=(8, 6), dpi_save=300)
#sc.pl.matrixplot(acts_epi, ctypes_dict, leiRes, dendrogram=True, standard_scale='var',
#                 colorbar_title='Z-scaled scores', cmap='RdBu_r',save=f'{Experiment}_{leiRes}_EpiORA_matrix.png',show=False)
#sc.set_figure_params(figsize=(6, 6), dpi_save=300)
#annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()
#adata_epi_NoSV.obs['TypeDef'] = [annotation_dict[clust] for clust in adata_epi_NoSV.obs[leiRes]]
#sc.pl.umap(adata_epi_NoSV,color=f'TypeDef',alpha=0.7,title=f'{Experiment}_Epithelial_Types',save=f'{Experiment}_{leiRes}_EpiORA_TypeDef_UMAP.png',show=False)






#FULL
markers_long_list = []
for gset_name in gset_df.columns:
    for gen in list(gset_df[gset_name].dropna()):
        markers_long_list.append([gen,gset_name])

markers_long = pd.DataFrame(markers_long_list)
col_names = ['genesymbol', 'cell_type']
markers_long.columns = col_names

dc.run_ora(
    mat=adata,
    net=markers_long,
    source='cell_type',
    target='genesymbol',
    min_n=3,
    verbose=True,
    use_raw=False
)

acts = dc.get_acts(adata, obsm_key='ora_estimate')
# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

for lineage in ['Epi','Str','Imm']:
    sc.pl.umap(acts, color=[col for col in gset_df.columns if lineage in col],save=f'{Experiment}_{lineage}_Full_ORA.png',size=300000/acts.n_obs, cmap='RdBu_r',show=False,vmin=0,vmax=12)

for lineage in ['Epi','Str','Imm']:
    sc.pl.umap(adata, color=[col for col in gset_df.columns if lineage in col],save=f'{Experiment}_{lineage}_Full_ORA.png',size=300000/acts.n_obs, cmap='RdBu_r',show=False,vmin=0,vmax=12)

#marker_dict = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']
#log_bool=False
#
#fig, axes = plt.subplots(1, len(marker_dict), figsize=(len(marker_dict)*5, 5))
#    # Plot histograms for each column
#for i, marker_name in enumerate(list(marker_dict)):
#    axes[i].hist(acts.obsm['ora_estimate'][marker_name], bins=100, color='skyblue',log=log_bool) #, ,range=(-0.2, 0.6)
#
#    # Add labels and title
#    axes[i].set_xlabel(marker_name + ' Score')
#    if log_bool:
#        axes[i].set_ylabel('Log Number of cells')
#    else:
#        axes[i].set_ylabel('Number of cells')
#    axes[i].set_title(marker_name + ' ORA Score Distribution')
#
#plt.tight_layout()
#plt.savefig(f'{Experiment}_ORAScore_Hist.pdf')
#plt.close()
leiRes='leiden_res3'
sc.set_figure_params(figsize=(8, 6), dpi_save=300)
for subtype in epi_marker_names_NoSV:
    acts.obs[subtype] = acts.obsm['ora_estimate'][subtype].copy()
    mean_values = acts.obs.groupby(leiRes,observed=False)[subtype].median()
    sorted_clusters = mean_values.sort_values(ascending=False).index.tolist()
    acts.obs.drop(columns=subtype, inplace=True)
    sc.pl.violin(acts, keys=subtype, groupby=leiRes,order=sorted_clusters,save=f'{Experiment}_{subtype}_{leiRes}_ORA_violin.png',show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

df = dc.rank_sources_groups(acts, groupby=leiRes, reference='rest', method='t-test_overestim_var')

n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
sc.set_figure_params(figsize=(8, 6), dpi_save=300)
sc.pl.matrixplot(acts, ctypes_dict, leiRes, dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r',save=f'{Experiment}_{leiRes}_ORA_matrix.png',show=False)
sc.set_figure_params(figsize=(6, 6), dpi_save=300)

annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()
adata.obs['TypeDef'] = [annotation_dict[clust] for clust in adata.obs[leiRes]]
sc.pl.umap(adata,color=f'TypeDef',alpha=0.7,title=f'{Experiment}_Epithelial_Types',save=f'{Experiment}_{leiRes}_FullORA_TypeDef_UMAP.png',show=False)

#sc.pl.umap(adata,size=size_list,color=f'TypeDef',alpha=0.7,title=f'{Experiment}_Epithelial_Types',save=f'{Experiment}_Epithelial_TypeDef_UMAP.png',show=False)





#Cluster rather than individual cells

#VENN 3 SUBTYPES
for gset_name in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1']:
    adata_epi.obs[f'{gset_name}_Type'] = [adata.obs[leiRes].isin(get_highscore_clusters(adata_epi,gset_name,leires=leiRes,medianthresh=0.65))][0]

#for gset_name in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1']:
#    adata_epi.obs[f'{gset_name}_Type'] = adata_epi.obs[gset_name] > threshold_genescoreSub

from matplotlib_venn import venn3
class_sets = []
for gset_name in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1']:
    class_sets.append(set(adata_epi.obs.index[adata_epi.obs[f'{gset_name}_Type']]))

venn_diagram = venn3(class_sets, set_labels=['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1'])
plt.title("Overlap in Epithelial Luminal Subtypes")
plt.tight_layout()
plt.savefig(f'{Experiment}_Epi_Subtypes_Clus_Venn3.pdf')
plt.close()


setL1 = set(adata_epi.obs.index[adata_epi.obs[f'Epi_Luminal_1_Type']])
setL2 = set(adata_epi.obs.index[adata_epi.obs[f'Epi_Luminal_2Psca_Type']])
setb = set(adata_epi.obs.index[adata_epi.obs[f'Epi_Basal_1_Type']])
setAmbig = set(setb & setL1 | setb & setL2 | setL1 & setL2)

dicTypeDef = {}
adata_epi.obs[f'TypeDef'] = pd.Categorical(['Ambiguous']*adata_epi.n_obs, categories=['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1','Ambiguous'])
for gset_name in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1']:
    dicTypeDef[gset_name] = set(adata_epi.obs.index[adata_epi.obs[f'{gset_name}_Type']]) - setAmbig
    adata_epi.obs.loc[adata_epi.obs.index.isin(dicTypeDef[gset_name]), 'TypeDef'] = gset_name

vminEpi = np.min([adata_epi.obs[gse] for gse in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1']])
vmaxEpi = np.max([adata_epi.obs[gse] for gse in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1']])
for gset_name in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Basal_1']:
    sc.pl.umap(adata_epi,color=gset_name,vmin=vminEpi,vmax=vmaxEpi,title=f'{Experiment}_{gset_name}_Score',save=f'{Experiment}_{gset_name}_UMAP.png',show=False)
sc.pl.umap(adata_epi,color=f'TypeDef',palette=['#a6cee3','#b2df8a','#fb9a99','#8a6be2','#bdbdbd'],title=f'{Experiment}_Epithelial_Type',save=f'{Experiment}_Epithelial_TypeDef_UMAP.png',show=False)


sc.pl.umap(adata,color=leiRes,legend_loc="on data",save=f"{Experiment}_{leiRes}",title=Experiment,show=False)
sc.pl.umap(adata_epi,color=leiRes,legend_loc="on data",save=f"{Experiment}_EpiOnly_{leiRes}",title=Experiment,show=False)

boxplot_gene_sigs(adata_epi,'Epi_Luminal_1','RNA_Epi_L1',clusters=leiRes)
boxplot_gene_sigs(adata_epi,'Epi_Luminal_2Psca','RNA_Epi_L2',clusters=leiRes)
boxplot_gene_sigs(adata_epi,'Epi_Basal_1','RNA_Epi_Bas',clusters=leiRes)


#Redoing clustering: maybe?
sc.pp.highly_variable_genes(adata_epi, layer="scran_normalization")
adata_epi.X = adata_epi.layers["scran_normalization"]
sc.pp.pca(adata_epi, svd_solver="arpack", use_highly_variable=True)
sc.pp.neighbors(adata_epi, use_rep="X_pca",n_neighbors=30, n_pcs=20)
sc.tl.umap(adata_epi)

for resLeiden in [.1,.15,.2]:
    print(f'Leiden clustering at {resLeiden} resolution')
    sc.tl.leiden(adata_epi, key_added=f"leiden_res{resLeiden}", resolution=resLeiden)
    sc.pl.umap(adata_epi,color=f'leiden_res{resLeiden}',legend_loc="on data",save=f"{Experiment}_EpiOnlyNewCluster_leiden_res{resLeiden}",title=Experiment,show=False)

boxplot_gene_sigs(adata_epi,'Epi_Luminal_1','RNA_Epi_L1',clusters='leiden_res0.1')
boxplot_gene_sigs(adata_epi,'Epi_Luminal_2Psca','RNA_Epi_L2',clusters='leiden_res0.1')
boxplot_gene_sigs(adata_epi,'Epi_Basal_1','RNA_Epi_Bas',clusters='leiden_res0.1')
sc.pl.umap(adata_epi,color=f'TypeDef',palette=['#a6cee3','#b2df8a','#fb9a99','#8a6be2','#bdbdbd'],title=f'{Experiment}_Epithelial_Type',save=f'{Experiment}_EpithelialNewClsuter_TypeDef_UMAP.png',show=False)
for gset_name in ['Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Luminal_3Foxi1','Epi_Basal_1']:
    sc.pl.umap(adata_epi,color=gset_name,vmin=vminEpi,vmax=vmaxEpi,title=f'{Experiment}_{gset_name}_Score',save=f'{Experiment}_{gset_name}NewCluster_UMAP.png',show=False)

adata_epi.write(f'{Experiment}_Epi{annot_ext}')





def Subtypes_Thresh_ScatterL3(adat,name,Experiment=Experiment,thresholds = [i / 100 for i in range(0, 151)]):
# distribution of uniquely called subtypes depending on genescoring thresholds
    l1 = []
    l2 = []
    l3 = []
    b = []
    ambig = []
    for th in thresholds:
        setL1 = set(adat.obs.index[adat.obs['Epi_Luminal_1'] > th])
        setL2 = set(adat.obs.index[adat.obs['Epi_Luminal_2Psca'] > th])
        setL3 = set(adat.obs.index[adat.obs['Epi_Luminal_3Foxi1'] > th])
        setb = set(adat.obs.index[adat.obs['Epi_Basal_1'] > th])
        setAmbig = set(setb & setL1 | setb & setL2 | setb & setL3 | setL1 & setL2 | setL1 & setL3 | setL2 & setL3)

        l1.append(len(setL1 - (setL2 | setb | setL3)))
        l2.append(len(setL2 - (setL1 | setb | setL3)))
        l3.append(len(setL3 - (setL1 | setb | setL2)))
        b.append(len(setb - (setL1 | setL2 | setL3)))
        ambig.append(len(setAmbig))

    plt.scatter(x=thresholds,y=l1, alpha=0.8, color='#a6cee3', edgecolor='black', linewidth=0.2, label="L1 only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=l2, alpha=0.8, color='#b2df8a', edgecolor='black', linewidth=0.2, label="L2 only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=l3, alpha=0.8, color='#fb9a99', edgecolor='black', linewidth=0.2, label="L3 only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=b, alpha=0.8, color='#8a6be2', edgecolor='black', linewidth=0.2, label="Basal only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=ambig, alpha=0.8, color='#bdbdbd', edgecolor='black', linewidth=0.2, label="Ambig",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.yscale('log')
    plt.ylabel('Log Number of Epithelial cells')
    plt.xlabel('Gene score threshold for subtype assignment')
    plt.title('Uniquely assigned Epithelial subtypes')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{Experiment}_{name}_NCells_UniqueSubtypes_Scatter.pdf')
    plt.close()

# distribution of uniquely called subtypes depending on genescoring thresholds
thresholds = [i / 100 for i in range(0, 151)]
#if 'ATAC' in DirRNA:
#    thresholds = [i / 100 for i in range(0, 40)]
l1 = []
l2 = []
l3 = []
b = []
ambig = []

reslist = []
for i,th in enumerate(thresholds):
    setL1 = set(get_highscore_clusters(adata_epi,'Epi_Luminal_1',leires='leiden_res3',medianthresh=th))
    setL2 = set(get_highscore_clusters(adata_epi,'Epi_Luminal_2Psca',leires='leiden_res3',medianthresh=th))
    setL3 = set(get_highscore_clusters(adata_epi,'Epi_Luminal_3Foxi1',leires='leiden_res3',medianthresh=th))
    setb = set(get_highscore_clusters(adata_epi,'Epi_Basal_1',leires='leiden_res3',medianthresh=th))
    setAmbig = set(setb & setL1 | setb & setL2 | setb & setL3 | setL1 & setL2 | setL1 & setL3 | setL2 & setL3)

    l1.append(len(setL1 - (setL2 | setb | setL3)))
    l2.append(len(setL2 - (setL1 | setb | setL3)))
    l3.append(len(setL3 - (setL1 | setb | setL2)))
    b.append(len(setb - (setL1 | setL2 | setL3)))
    ambig.append(len(setAmbig))

    reslist.append((th,l1[i]+l2[i]+l3[i]))

plt.scatter(x=thresholds,y=l1, alpha=0.8, color='#a6cee3', edgecolor='black', linewidth=0.2, label="L1 only",s=plt.rcParams['lines.markersize'] ** 1.5)
plt.scatter(x=thresholds,y=l2, alpha=0.8, color='#b2df8a', edgecolor='black', linewidth=0.2, label="L2 only",s=plt.rcParams['lines.markersize'] ** 1.5)
plt.scatter(x=thresholds,y=l3, alpha=0.8, color='#fb9a99', edgecolor='black', linewidth=0.2, label="L3 only",s=plt.rcParams['lines.markersize'] ** 1.5)
plt.scatter(x=thresholds,y=b, alpha=0.8, color='#8a6be2', edgecolor='black', linewidth=0.2, label="Basal only",s=plt.rcParams['lines.markersize'] ** 1.5)
plt.scatter(x=thresholds,y=ambig, alpha=0.8, color='#bdbdbd', edgecolor='black', linewidth=0.2, label="Ambig",s=plt.rcParams['lines.markersize'] ** 1.5)
plt.ylabel(f'Number of {leiRes} leiden cluster')
plt.xlabel('Gene score threshold for subtype attribution')
plt.title('Uniquely assigned Epithelial subtypes')
plt.legend()
plt.tight_layout()
plt.savefig(f'{Experiment}_Epi_NClusters_in_Subtypes_Thresh_Scatter.pdf')
plt.close()


# single gene Scores? lol
#adata.var.index = adata.var.index.str.upper()
#main_markers = {"Epithelial": ["Epcam"], "Immune": ["Ptprc"]}  # present
#for gset_name in main_markers.keys():
#    gs = main_markers[gset_name]
#    sc.tl.score_genes(adata, gs, ctrl_size=25, n_bins=25, score_name=gset_name,use_raw=False)
#    sc.pl.umap(adata,color=gset_name,save=f'{Experiment}_{gset_name}_AnnotMain.png',title=f'{Experiment}_{gset_name}',show=False)
#hist_gene_sigs(adata,main_markers,'Main',threshold_genescoreMain)


## threshold based on the Distribution
#adata.obs['Epithelial_Type'] = adata.obs['Epithelial'] > threshold_genescoreMain
#adata.obs['Immune_Type'] = adata.obs['Immune'] > threshold_genescoreMain
#
##VENN IMM vs EPI 
#from matplotlib_venn import venn2
## Count occurrences of each category combination
#epithelial_only = adata.obs[(adata.obs['Epithelial_Type']) & (~adata.obs['Immune_Type'])]
#immune_only = adata.obs[(~adata.obs['Epithelial_Type']) & (adata.obs['Immune_Type'])]
#both = adata.obs[(adata.obs['Epithelial_Type']) & (adata.obs['Immune_Type'])]
#
## Create the Venn diagram
#venn2(subsets=(len(epithelial_only), len(immune_only), len(both)),
#      set_labels=('Epithelial', 'Immune'))
#
## Display the plot
#plt.suptitle("Overlap between Epithelial and Immune Types")
#plt.title(f"Total cells: {adata.n_obs}")
#plt.tight_layout()
#plt.savefig(f'{Experiment}_Epi_Imm_Venn2.pdf')
#plt.close()
#
##Subset adata
#adata_epi = adata[(adata.obs['Epithelial_Type']) & ~(adata.obs['Epithelial_Type'] & adata.obs['Immune_Type'])].copy()
##adata_imm = adata[(adata.obs['Immune_Type']) & ~(adata.obs['Epithelial_Type'] & adata.obs['Immune_Type'])].copy()






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

plt.title("RNA Gene Scores")
plt.tight_layout()
plt.savefig(Experiment+'_RNAGeneScore_Hist.pdf')
plt.close()

#UMAPS
for gset_name in gset_dict.keys():
    print(f'{gset_name}: {len(gs)}genes Karthaus2020 TableS8')
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


#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Imm*.png --outfile SING_Imm.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Str*.png --outfile SING_Str.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Epi*.png --outfile SING_Epi.pdf
#
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Imm*.png --outfile ENZ_Imm.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Str*.png --outfile ENZ_Str.pdf
#pdfjam --nup 4x2 --landscape --a4paper figures/umap*Epi*.png --outfile ENZ_Epi.pdf




















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


#Celltyping based on Paper genes 
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L1']], ctrl_size=50, n_bins=25, score_name='L1_genes_Paper')
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L2']], ctrl_size=50, n_bins=25, score_name='L2_genes_Paper')
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L3']], ctrl_size=50, n_bins=25, score_name='L3_genes_Paper')
sc.tl.score_genes(adata, [m.upper() for m in markerGenes['L']], ctrl_size=50, n_bins=25, score_name='L_genes_Paper')

sc.pl.umap(adata,color=["L1_genes_Paper"],save=Experiment+'_L1_genes_fromMain',title=f'{Experiment}_L1 sig paper',show=False)
sc.pl.umap(adata,color=["L2_genes_Paper"],save=Experiment+'_L2_genes_fromMain',title=f'{Experiment}_L2 sig paper',show=False)
sc.pl.umap(adata,color=["L3_genes_Paper"],save=Experiment+'_L3_genes_fromMain',title=f'{Experiment}_L3 sig paper',show=False)
sc.pl.umap(adata,color=["L_genes_Paper"],save=Experiment+'_L_genes_fromMain',title=f'{Experiment}_L sig paper',show=False)


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


