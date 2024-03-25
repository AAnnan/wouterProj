import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
import anndata as ad
import os
sc.set_figure_params(figsize=(6, 6), dpi_save=300)
print("Scanpy version:",sc.__version__)

### TO CHECK WHEN CHANGING SAMPLES ###
Experiment='Wouter21_ENZ_CB'
### TO CHECK WHEN CHANGING SAMPLES ###

refDir = '/mnt/etemp/ahrmad/wouter/refs'
#refDir = './refs'

if not os.path.exists(refDir):
    raise ValueError('Check ref folder paths')
gene_set_table = f'{refDir}/table_s8_summary.txt'


post_ext = '_post.h5ad'
annot_ext = '_Annotated.h5ad'

leiRes = 'leiden_res1'


#adata import
adata = sc.read_h5ad(Experiment + post_ext)
sc.pl.umap(adata,color=leiRes,legend_loc="on data",save=f"{Experiment}_{leiRes}",title=Experiment,show=False)



### Top Level Lineages (Epi,Imm,Strm)

#Hist of the Distribution
def hist_SoloGenesScore(adata,genes,name_genelist,raw=False,Experiment=Experiment,log_bool=True):
    lay = 'scran_normalization'
    if raw:
        lay = 'No_normalization'
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
    plt.savefig(f'{Experiment}_{name_genelist}Score_Hist.pdf')
    plt.close()

hist_SoloGenesScore(adata,["Epcam","Ptprc"],'Epcam_Ptprc')
adata.obs['EpcamScran'] = adata.layers['scran_normalization'].T[adata.var.index == 'Epcam'].toarray()[0]
adata.obs['PtprcScran'] = adata.layers['scran_normalization'].T[adata.var.index == 'Ptprc'].toarray()[0]

sc.pl.umap(adata,color='EpcamScran',save=f'{Experiment}_EpcamScranScore.png',title=f'{Experiment}_Epcam',show=False,layer='scran_normalization', vmin=0, vmax=2)
sc.pl.umap(adata,color='PtprcScran',save=f'{Experiment}_PtprcScranScore.png',title=f'{Experiment}_Ptprc',show=False,layer='scran_normalization', vmin=0, vmax=2)
adata.obs['EpcamExpressing'] = pd.Categorical(adata.layers['scran_normalization'].T[adata.var.index == 'Epcam'].toarray()[0] > 0, categories=[True,False])
sc.pl.umap(adata,color='EpcamExpressing',save=f'{Experiment}_EpcamExpressing.png',title=f'{Experiment}_EpcamExpressing',show=False)


#Hist of the Distribution
def hist_gene_sigs(adata,marker_dict,name_dict,thresh_gs,Experiment=Experiment,log_bool=True):
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
    plt.savefig(f'{Experiment}_{name_dict}Score_Hist.pdf')
    plt.close()
    return 0

def boxplot_gene_sigs(mod,gsetName,modality,Experiment=Experiment,clusters='leiden_res0.25'):
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
    plt.savefig(f'{Experiment}_{modality}boxplots_{clusters}_{gsetName}Score.pdf')
    plt.close()

boxplot_gene_sigs(adata,'EpcamScran','RNA',clusters=leiRes)
boxplot_gene_sigs(adata,'PtprcScran','RNA',clusters=leiRes)


def get_highscore_clusters(adata_obj,gene_obs,leires='leiden_res2.5',medianthresh=0):
    HighScoreClusters = []
    for clus in adata_obj.obs[leires].cat.categories:
        MedClus = np.median(adata_obj.obs[gene_obs][adata_obj.obs[leires]==clus])
        if MedClus > medianthresh:
            HighScoreClusters.append(clus)
    return HighScoreClusters


HighEpcamClusters = get_highscore_clusters(adata,'EpcamScran',leires=leiRes,medianthresh=0)
HighPtprcClusters = get_highscore_clusters(adata,'PtprcScran',leires=leiRes,medianthresh=0)

assert len(np.intersect1d(HighEpcamClusters,HighPtprcClusters))==0,'Some clusters overlap between Epcam and Ptprc'

adata_epi = adata[adata.obs[leiRes].isin(HighEpcamClusters)].copy()
adata_imm = adata[adata.obs[leiRes].isin(HighPtprcClusters)].copy()
adata_str = adata[~adata.obs[leiRes].isin(HighEpcamClusters+HighPtprcClusters)].copy()

adata.obs['TopLevelLineages'] = 'Stromal'
adata.obs.loc[[adata.obs[leiRes].isin(HighPtprcClusters)][0],'TopLevelLineages'] = 'Immune'
adata.obs.loc[[adata.obs[leiRes].isin(HighEpcamClusters)][0],'TopLevelLineages'] = 'Epithelial'

sc.pl.umap(adata,color='TopLevelLineages',save=f'{Experiment}_{leiRes}_TopLevelLineages.png',title=f'{Experiment}_{leiRes}_Top_Level_Lineages',show=False)

#HIST cells in samples, epi vs total 
plt.hist(adata.obs['batch'], bins=len(adata.obs['batch'].unique()), alpha=0.5, label='All')
plt.hist(adata_epi.obs['batch'], bins=len(adata_epi.obs['batch'].unique()), alpha=0.5, label='Epithelial')
plt.xlabel(None)
plt.ylabel('Number of cells')
plt.title(f'Epithelial cells in {Experiment} Samples')
plt.legend()
plt.xticks(rotation=45, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.savefig(f'{Experiment}_{leiRes}_ProportionEpi_Barplot.pdf')
plt.close()





### Finer Cell Typing (Basal, L1, L2 esp)

#Epi subtypes Gene Scoring based on Supp gene lists
gset_df = pd.read_csv(gene_set_table, sep='\t')
epi_marker_names = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Luminal_3Foxi1','Epi_Basal_SV','Epi_Luminal_SV']
epi_gset_dict = dict()

for gset_name in epi_marker_names:
    epi_gset_dict[gset_name] = list(gset_df[gset_name].dropna())

# Clean up (removes 1-4 genes per module)
epi_gset_dict = {key: [item for item in value if item in adata_epi.var_names] for key, value in epi_gset_dict.items()}
#for key in epi_gset_dict:
#    print(f'{key}\n{len(epi_gset_dict2[key])}\n{len(epi_gset_dict[key])}')

# Gene scoring on epi subset
for gset_name in epi_gset_dict.keys():
    gs = epi_gset_dict[gset_name]
    print(f'{gset_name}: {len(gs)}genes Karthaus2020 TableS8')
    ctrl = len(gs)
    sc.tl.score_genes(adata_epi, gs, ctrl_size=ctrl, n_bins=25, score_name=gset_name,use_raw=False)

hist_gene_sigs(adata_epi,epi_marker_names,'Epithelial_Markers',thresh_gs=None)

allEpiscores = np.array([adata_epi.obs[gse].tolist() for gse in epi_marker_names])
vminEpi = np.min(allEpiscores)
vmaxEpi = np.max(allEpiscores)
for gset_name in epi_marker_names:
    sc.pl.umap(adata_epi,color=gset_name,vmin=vminEpi,vmax=vmaxEpi,title=f'{Experiment}_{gset_name}Score',save=f'{Experiment}_{gset_name}_EpiAnnotSuppS8.png',show=False)

SV_cells = set(adata_epi.obs.index[adata_epi.obs['Epi_Basal_SV'] >= threshold_SV]) | set(adata_epi.obs.index[adata_epi.obs['Epi_Luminal_SV'] >= threshold_SV])
adata_epi_NoSV = adata_epi[~adata_epi.obs.index.isin(SV_cells)].copy()

def Subtypes_Thresh_ScatterL2(adat,name,Experiment=Experiment,thresholds = [i / 100 for i in range(0, 151)]):
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
        setAmbig = set(setb & setL1 | setb & setL2 | setL1 & setL2)

        l1.append(len(setL1 - (setL2 | setb)))
        l2.append(len(setL2 - (setL1 | setb)))
        b.append(len(setb - (setL1 | setL2)))
        ambig.append(len(setAmbig))
        NoTypes.append(len(set(adat.obs.index) - set(setL1 | setL2 | setb)))

    plt.scatter(x=thresholds,y=np.array(l1)/n_cells, alpha=0.8, color='#a6cee3', edgecolor='black', linewidth=0.2, label="L1 only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(l2)/n_cells, alpha=0.8, color='#b2df8a', edgecolor='black', linewidth=0.2, label="L2 only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(b)/n_cells, alpha=0.8, color='#8a6be2', edgecolor='black', linewidth=0.2, label="Basal only",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(ambig)/n_cells, alpha=0.8, color='#bdbdbd', edgecolor='black', linewidth=0.2, label="Ambig",s=plt.rcParams['lines.markersize'] ** 1.5)
    plt.scatter(x=thresholds,y=np.array(NoTypes)/n_cells, alpha=0.8, color='#fb9a99', edgecolor='black', linewidth=0.2, label="Untyped",s=plt.rcParams['lines.markersize'] ** 1.5)
    #plt.yscale('log')
    plt.ylabel(f'% of  Epithelial cells (w/o SV)')
    plt.xlabel('Gene score threshold for subtype assignment')
    plt.title('Uniquely assigned Epithelial subtypes')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{Experiment}_{name}_NCells_UniqueSubtypes_Scatter.pdf')
    plt.close()

Subtypes_Thresh_ScatterL2(adata_epi_NoSV,'EPI_NoSV')


epi_marker_names_NoSV = ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']


if 'ENZ' in Experiment:
    threshold_SV = 0.70
    threshold_genescoreSub = 0.40
elif 'SING' in Experiment:
    threshold_SV = 0.70
    threshold_genescoreSub = 0.55

#VENN 3 SUBTYPES
for gset_name in epi_marker_names_NoSV:
    adata_epi_NoSV.obs[f'{gset_name}_Type'] = adata_epi_NoSV.obs[gset_name] > threshold_genescoreSub

class_sets = []
for gset_name in epi_marker_names_NoSV:
    class_sets.append(set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs[f'{gset_name}_Type']]))

venn_diagram = venn3(class_sets, set_labels=epi_marker_names_NoSV)
plt.title("Overlap in Epithelial Luminal Subtypes")
plt.tight_layout()
plt.savefig(f'{Experiment}_Epi_Subtypes_Venn3.pdf')
plt.close()


setL1 = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Luminal_1'] > threshold_genescoreSub])
setL2 = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Luminal_2Psca'] > threshold_genescoreSub])
setb = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs['Epi_Basal_1'] > threshold_genescoreSub])

All_Ambig = set(setb & setL1 | setb & setL2 | setL1 & setL2)
L1L2_Ambig = set(setL1 & setL2)
L1B_Ambig = set(setL1 & setb)
L2B_Ambig = set(setb & setL2)
L1L2B_Ambig = set(setb & setL1 & setL2)
setL1_NoAmbig = setL1 - All_Ambig
setL2_NoAmbig = setL2 - All_Ambig
setb_NoAmbig = setb - All_Ambig

adata_epi_NoSV.obs[f'TypeDef'] = pd.Categorical(['LowScores']*adata_epi_NoSV.n_obs, categories=epi_marker_names_NoSV+['L1/L2_Ambig','L1/B_Ambig','L2/B_Ambig','L1/L2/B_Ambig','LowScores'])
for gset_name,TypeDef_cells_gset_name in zip(epi_marker_names_NoSV+['L1/L2_Ambig','L1/B_Ambig','L2/B_Ambig','L1/L2/B_Ambig'],[setb_NoAmbig,setL1_NoAmbig,setL2_NoAmbig,L1L2_Ambig,L1B_Ambig,L2B_Ambig,L1L2B_Ambig]):
    adata_epi_NoSV.obs.loc[adata_epi_NoSV.obs.index.isin(TypeDef_cells_gset_name), 'TypeDef'] = gset_name

#TypeDef_cells_gset_name = set(adata_epi_NoSV.obs.index[adata_epi_NoSV.obs[gset_name] > threshold_genescoreSub]) - setAmbig
size_list = []
for cel in adata_epi_NoSV.obs.index:
    if adata_epi_NoSV.obs['TypeDef'][cel] in epi_marker_names_NoSV:
        size_list.append(360000/adata_epi_NoSV.n_obs)
    else:
        size_list.append(220000/adata_epi_NoSV.n_obs)

#sc.pl.umap(adata_epi_NoSV,color=f'TypeDef',title=f'{Experiment}_Epithelial_Types',save=f'{Experiment}_Epithelial_TypeDef_UMAP.png',show=False)
sc.pl.umap(adata_epi_NoSV,size=size_list,color=f'TypeDef',alpha=0.7,title=f'{Experiment}_Epithelial_Types',save=f'{Experiment}_Epithelial_TypeDef_UMAP.png',show=False)










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


