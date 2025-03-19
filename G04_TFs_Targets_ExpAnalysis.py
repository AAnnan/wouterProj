import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import matplotlib.pyplot as plt
from collections import Counter
from gseapy import Biomart
import os 

from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42

#from G03_Annotation_Functions import *

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
print("Scanpy version:",sc.__version__)

### TO CHECK WHEN CHANGING SAMPLES ###
Experiment='Wouter21_SING_CB'
### TO CHECK WHEN CHANGING SAMPLES ###
refDir = '/mnt/etemp/ahrmad/wouter/refs'
adata_dir = '/mnt/etemp/ahrmad/wouter/RNA_CB/annot_other_norm/'

if not os.path.exists(refDir):
    raise ValueError('Check ref folder paths')
gene_set_table = f'{refDir}/table_s8_summary.txt'
cell_cycle_gene_table = f'{refDir}/Mouse_G2M_S_genes_Tirosh2016.csv'

#Epi subtypes Supp gene lists
gset_df = pd.read_csv(gene_set_table, sep='\t')
cc_df = pd.read_csv(cell_cycle_gene_table)

post_ext = '_post.h5ad'
lum_ext = '_annot_Lum_newSigs.h5ad'
all_ext = '_annot_All_newSigs.h5ad'


#TF list
log2fc = 0.5
L1_CD28_RD3 = pd.read_csv(f'/mnt/etemp/ahrmad/wouter/refs/Wouter21_SING_TF_motifsHOCOMOCOv12_L1_TP_Day28_RegenDay3_minFC{log2fc}.csv')
L1_I_CD28 = pd.read_csv(f'/mnt/etemp/ahrmad/wouter/refs/Wouter21_SING_TF_motifsHOCOMOCOv12_L1_TP_Intact_Day28_minFC{log2fc}.csv')
L2_CD28_RD3 = pd.read_csv(f'/mnt/etemp/ahrmad/wouter/refs/Wouter21_SING_TF_motifsHOCOMOCOv12_L2_Ambig_L1_L2_TP_Day28_RegenDay3_minFC{log2fc}.csv')
L2_I_CD28 = pd.read_csv(f'/mnt/etemp/ahrmad/wouter/refs/Wouter21_SING_TF_motifsHOCOMOCOv12_L2_Ambig_L1_L2_TP_Intact_Day28_minFC{log2fc}.csv')

#NEW vv
#L1_I_CD28 = pd.read_csv(f'/mnt/etemp/ahrmad/wouter/Wouter21_SING/annot/figL1/figures/Wouter21_SING_TF_motifsHOCOMOCOv12_L1_TP_Intact_Day28_minFC1.tsv',sep='\t')
L1_I_CD28 = pd.read_csv(f'/mnt/etemp/ahrmad/wouter/Wouter21_SING/annot/20241214_L1/Motifs_Delta_H12.tsv',sep='\t')
num_ex=40
# Sort by absolute delta rank for each metric
df_log2_fc = L1_I_CD28.sort_values(by="delta_log2_fc", ascending=False).reset_index()
df_adj_pval = L1_I_CD28.sort_values(by="delta_rank_adj_pval", ascending=False).reset_index()

extreme_log2_fc = pd.concat([df_log2_fc.iloc[:num_ex]])
extreme_adj_pval = pd.concat([df_adj_pval.iloc[:num_ex]])
extreme_points=pd.concat([extreme_log2_fc, extreme_adj_pval])

TFs_dict_L1 = {
    'L1_I_CD28': extreme_points['id'].str.split('.').str[0].unique().tolist(),
}
#NEW ^^

TFs_dict_L1 = {
    'L1_CD28_RD3': L1_CD28_RD3['id'].str.split('.').str[0].unique().tolist(),
    'L1_I_CD28': L1_I_CD28['id'].str.split('.').str[0].unique().tolist(),
}

TFs_dict_L2 = {
    'L2_CD28_RD3': L2_CD28_RD3['id'].str.split('.').str[0].unique().tolist(),
    'L2_I_CD28': L2_I_CD28['id'].str.split('.').str[0].unique().tolist()
}

human_to_mouse_gene_map = {
    'ANDR': 'Ar',
    'ATF6A': 'Atf6',
    'EGR4': 'Egr4',
    'FEV': 'Fev',
    'GCR': 'Nr3c1',
    'KAISO': 'Zbtb33',
    'KLF14': 'Klf14',
    'LRRF1': 'Lrrfip1',
    'MCR': 'Mc2r',
    'MYCN': 'Mycn',
    'PRGR': 'Pgr',
    'SRBP1': 'Srebf1',
    'TYY1': 'Tcf3',
    'TYY2': 'Tcf4',
    'ZBED1': 'Zbed1',
    'ZBT14': 'Zbtb14',
    'ZBT37': 'Zbtb37',
    'ZN100': 'Zfp100',
    'ZN141': 'Zfp141',
    'ZN267': 'Zfp267',
    'ZN407': 'Zfp407',
    'ZN417': 'Zfp317',
    'ZN441': 'Zfp441',
    'ZN454': 'Zfp454',
    'ZN460': 'Zfp460',
    'ZN519': 'Zfp819',
    'ZN543': 'Zfp543',
    'ZN610': 'Zfp610',
    'ZN611': 'Zfp611',
    'ZN704': 'Zfp704',
    "BHE41": "Dec2",
    "NGN1": "Neurog1",
    "PKNX1": "Pknox1",
    "UNC4": "Uncx",
    "ZN649": "Zfp649",
    "BHE40": "BHLHE40",
    "CR3L4": "Creb3l4",
    "DPRX": "Dprx",
    "HXD12": "Hoxd-12",
    "MGAP": "Rapgef1",
    "SRBP2": "Srebf2",
    "ZBT17": "Zbtb17",
    "ZBT7B": "Zbtb7b",
    "ZBT7C": "Zbtb7c",
    "ZN211": "Zfp211",  # Zinc finger protein 211
    "NGN2": "Neurog2",  # Neurogenin 2
    "PO5F1": "Pou5f1",  # POU domain, class 5, transcription factor 1
    "BARH2": "Barhl2",  # BarH-like homeobox 2
    "P73": "Trp73",     # Transformation-related protein 73
    "ZN621": "Zfp621",  # Zinc finger protein 621
    "ZN565": "Zfp565",  # Zinc finger protein 565
}

TopHits_250219 = [
    "Zbtb33",  # ZBTB33 (Kaiso)
    "Creb3l4",  # CREB3L4
    "Nfya",     # NFYA
    "Creb1",    # CREB1
    "Atf3",     # ATF3
    "Jdp2",     # JDP2
    "Crem",     # CREM
    "Nfyb",     # NFYB
    "Atf7",     # ATF7
    "Nfyc",     # NFYC
    "Nrf1",     # NRF1
    "Tfe3",     # TFE3
    "Usf2",     # USF2
    "Arnt2",    # ARNT2
    "Mlxipl",   # MLXIPL (ChREBP)
    "Yy1",      # YY1
    "Yy2",      # YY2
    "Nr3c2",    # NR3C2 (Mineralocorticoid receptor)
    "Nr3c1",    # NR3C1 (Glucocorticoid receptor)
    "Ar",       # AR (Androgen receptor)
    "Pgr",      # PGR (Progesterone receptor)
    "Gata1",    # GATA1
    "Gata6",    # GATA6
    "Gata2",    # GATA2
    "Trps1",    # TRPS1
    "Gata4",    # GATA4
    "Gata3",    # GATA3
    "Gata5",    # GATA5
    "Foxa1",    # FOXA1
    "Cdx2",     # CDX2
    "Dmrt1",    # DMRT1
    "Tal1",     # TAL1
    "Sox21",    # SOX21
    "Foxl1",    # FOXL1
    "Hmga1"     # HMGA1
]


adata = sc.read_h5ad(adata_dir + Experiment + lum_ext)

for gene in TopHits_250219:
    print(gene, gene in adata.var_names)
#Missing: Sox21,Cdx2,Gata5,Yy2

# ADD FAKE GENES (EXP 0)
# Genes to replace and their replacements
genes_to_replace = ['Yod1', 'Yme1l1', 'Ypel4', 'Ywhab']
new_genes = ['Sox21', 'Cdx2', 'Gata5', 'Yy2']
rename_dict = dict(zip(genes_to_replace, new_genes))

# Check if the genes to replace are in adata.var_names
for gene in genes_to_replace:
    if gene in adata.var_names:
        # Get the index of the gene
        gene_index = adata.var_names.get_loc(gene)
        
        # Set the expression of the gene to zero
        adata.X[:, gene_index] = 0

# Dictionary mapping old names to new names

# Rename genes
adata.var_names = adata.var_names.to_series().replace(rename_dict)

# Ensure var_names are unique
adata.var_names_make_unique()
# ADD FAKE GENES (EXP 0) ^

############################
#### Functions
############################

def get_MM_homologs_onlyMM(aliases,adat_f,human_to_mouse_gene_map=human_to_mouse_gene_map):
    if not aliases:
        print(f'Number of Human genes selected for conversion: {len(aliases)}')
        return []

    bm = Biomart()
    queries ={'external_gene_name': aliases} # need to be a dict object
    results = bm.query(dataset='hsapiens_gene_ensembl',
                       attributes=['ensembl_gene_id','external_gene_name',
                           'mmusculus_homolog_ensembl_gene',
                           'mmusculus_homolog_associated_gene_name'],
                       filters=queries)

    if results is None:
        print(f'Number of Human genes selected for conversion: {len(aliases)}')
        print(f'Number of Human genes not found in DB: {len(aliases)}')
        print(f'Human genes not found in DB \n{aliases}')
        return []

    results = results.dropna(subset=['mmusculus_homolog_associated_gene_name'])

    mm_genes = results['mmusculus_homolog_associated_gene_name']
    mm_genes_not_in_adata = mm_genes[~mm_genes.isin(adat_f.var.index)].tolist()
    mm_genes_in_adata = mm_genes[mm_genes.isin(adat_f.var.index)].tolist()
    hs_genes_not_in_biomart = [g for g in aliases if g not in results['external_gene_name'].tolist()]

    print(f'Number of Human genes selected for automatic conversion: {len(aliases)}')
    print(f'Number of Mouse homologs not found in ADATA: {len(mm_genes_not_in_adata)}/{len(aliases)}')
    if len(mm_genes_not_in_adata) > 0:
        print(f'Mouse homologs not found in ADATA (removed):\n{mm_genes_not_in_adata}')

    print(f'Number of Human genes not found in DB: {len(hs_genes_not_in_biomart)}/{len(aliases) - len(mm_genes_not_in_adata)}')
    if len(hs_genes_not_in_biomart) > 0:
        print(f'Human genes not found in DB (will do manual conversion):\n{hs_genes_not_in_biomart}')

    # Initialize list for mouse gene symbols
    add_mouse_genes = []
    # Track human genes that are not found in the mapping
    not_found_genes = []

    # Loop through each human gene and map it to the corresponding mouse gene
    for gene in hs_genes_not_in_biomart:
        if gene in human_to_mouse_gene_map:
            add_mouse_genes.append(human_to_mouse_gene_map[gene])
        else:
            not_found_genes.append(gene)
    
    # If there are missing genes, output a warning
    if not_found_genes:
        print(f"Warning: The following human genes were not found in the conversion list: {', '.join(not_found_genes)}")
    
    add_mouse_genes = pd.Series(add_mouse_genes)
    add_mouse_genes_in_adata = add_mouse_genes[add_mouse_genes.isin(adat_f.var.index)].tolist()
    add_mouse_genes_not_in_adata = add_mouse_genes[~add_mouse_genes.isin(adat_f.var.index)].tolist()
    
    print("Manually Converted mouse genes in ADATA:\n", add_mouse_genes_in_adata)
    if len(add_mouse_genes_not_in_adata) > 0:
        print("Manually Converted mouse genes not in ADATA (removed):\n", add_mouse_genes_not_in_adata)

    return mm_genes_in_adata + add_mouse_genes_in_adata

def get_MM_homologs(aliases,adat_f):
    bm = Biomart()
    queries ={'external_gene_name': aliases} # need to be a dict object
    results = bm.query(dataset='hsapiens_gene_ensembl',
                       attributes=['ensembl_gene_id','external_gene_name',
                           'mmusculus_homolog_ensembl_gene',
                           'mmusculus_homolog_associated_gene_name'],
                       filters=queries)
    print('before rem NA ',len(results))
    results = results.dropna(subset=['mmusculus_homolog_associated_gene_name'])
    print('After rem NA ',len(results))

    results = results[results['mmusculus_homolog_associated_gene_name'].isin(adat_f.var.index)]
    print('After check adata ',len(results))
    return results['mmusculus_homolog_associated_gene_name'].tolist(),(results['external_gene_name'] + ' / ' + results['mmusculus_homolog_associated_gene_name']).tolist()


def annot_DotPlot(adat_f,TF_ana_f,cat_order_f,gene_f,Experiment=Experiment):
        typ_list = adat_f.obs['Annotation'].cat.categories.tolist()
        # Extract the expression data for 'Atf1' from adata_inUse (can be .X or specific layer)
        gene_expr = adat_f[:, gene_f].X.copy()

        stack_list = []
        for typeCell in typ_list:
            mask_typ = np.array(adat_f.obs['Annotation'] == typeCell)[:, np.newaxis]

            # Create new arrays for each type with non-matching cells set to 0
            gene_var = np.where(mask_typ, gene_expr, 0)

            stack_list.append(gene_var)

        # Concatenate the gene-specific arrays into one matrix
        new_expr_data = np.hstack(stack_list)

        # Create a new AnnData object
        adata_new = ad.AnnData(X=new_expr_data, 
            obs=pd.DataFrame({'timePoint':adat_f.obs['timePoint'].copy()}), var={'new_names': typ_list})

        # Now adata_new contains the modified gene expression profiles

        sc.pl.DotPlot.DEFAULT_SMALLEST_DOT = 0.2
        sc.pl.DotPlot.DEFAULT_DOT_EDGELW = 0.12
        sc.pl.DotPlot.DEFAULT_SIZE_LEGEND_TITLE = f"Fraction of\n{gene_f} exp. cells (%)"
        sc.pl.DotPlot.DEFAULT_COLOR_LEGEND_TITLE = f"Mean expression"
        sc.pl.DotPlot.DEFAULT_LEGENDS_WIDTH = 2.2  # inches
        sc.pl.DotPlot.DEFAULT_PLOT_X_PADDING = 0.6  # a unit is the distance between two x-axis ticks
        sc.pl.DotPlot.DEFAULT_PLOT_Y_PADDING = 0.6  # a unit is the distance between two y-axis ticks
        
        suf = 'Abs_AllTimePoints' if len(cat_order_f)>3 else 'Abs'
        sc.pl.DotPlot(adata_new, var_names=typ_list, groupby='timePoint', use_raw=False,
            mean_only_expressed=True, expression_cutoff=0,standard_scale=None,#'var',  #defaults stuff
            categories_order=cat_order_f,gene_symbols='new_names',cmap='viridis',
            title=f'{Experiment}').add_totals(color='lightpink').swap_axes().savefig(f'./figures/{Experiment}_{TF_ana_f}_{gene_f}_{suf}.pdf')
        plt.close()

        return 0

def fast_dotplot(adat_f,TF_targets_f,TF_ana_f,cat_order_f,suf_f,Experiment=Experiment):

        if 'Rel' in suf_f:
            std_scale = 'var'
        elif 'Abs' in suf_f:
            std_scale = None

        if 'Top' in TF_ana_f:
            leg1 = TF_ana_f
        else:
            info = TF_ana_f.split('_')
            leg1 = info[0]
            leg2=f'\n {info[1]} vs {info[2]} DiffAcc TFs'
        sc.pl.DotPlot.DEFAULT_LEGENDS_WIDTH = 2.2  # inches

        sc.pl.DotPlot(adat_f, var_names=TF_targets_f, groupby='timePoint', use_raw=False,
            mean_only_expressed=True, expression_cutoff=0,standard_scale=std_scale,#'var', #defaults stuff
            categories_order=cat_order_f,cmap='bwr',
            title=f'{Experiment} {leg1}').add_totals(color='lightpink').savefig(f'./figures/{Experiment}_{TF_ana_f}_{suf_f}.pdf')
        plt.close()
        return 0


# Get Mouse 
for dictspl in [TFs_dict_L1]:
    for k in dictspl.keys():
        print(f'\n\n{k}')
        dictspl[k] = get_MM_homologs_onlyMM(dictspl[k],adata)

# Subset only per Luminal Types
cat_order_all = pd.Categorical([ 'Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'], categories=['Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'])
cat_order_top = pd.Categorical([ 'Intact','Day28', 'RegenDay3'], categories=['Intact','Day28', 'RegenDay3'])
cat_order_top = pd.Categorical([ 'Intact','Day28'], categories=['Intact','Day28'])

    # In case we want to retrieve back the human aliases
    #TF_targets,dot_names = get_MM_homologs(TFs_dict_L1[TF_ana],adata_inUse)
    #name_mapping = dict(zip(TF_targets,dot_names))
    #adata_inUse.var['dot_names'] = adata_inUse.var.index.to_series().replace(name_mapping)

    #in sc dotplot:
    #put var_names=dot_names and gene_symbols="dot_names",


#L1
adata_L1 = adata[np.isin(adata.obs['Annotation'],['L1'])].copy()
adata_L1_topTP = adata_L1[np.isin(adata_L1.obs['timePoint'],cat_order_top.tolist())].copy()


TF_targets = [g for g in TopHits_250219 if g in adata_L1.var.index]
fast_dotplot(adata_L1,TF_targets,"TopHits_250219",cat_order_all,'Abs_AllTP')
fast_dotplot(adata_L1,TF_targets,"TopHits_250219",cat_order_all,'Rel_AllTP')

TF_targets = [g for g in TopHits_250219 if g in adata_L1_topTP.var.index]
fast_dotplot(adata_L1_topTP,TF_targets,"TopHits_250219",cat_order_top,'Abs_TopTP')    
fast_dotplot(adata_L1_topTP,TF_targets,"TopHits_250219",cat_order_top,'Rel_TopTP')


for TF_ana in TFs_dict_L1.keys():
    print(TF_ana)

    TF_targets = [g for g in TFs_dict_L1[TF_ana] if g in adata_L1.var.index]
    fast_dotplot(adata_L1,TF_targets,TF_ana,cat_order_all,'Abs_AllTP')
    fast_dotplot(adata_L1,TF_targets,TF_ana,cat_order_all,'Rel_AllTP')
    
    TF_targets = [g for g in TFs_dict_L1[TF_ana] if g in adata_L1_topTP.var.index]
    fast_dotplot(adata_L1_topTP,TF_targets,TF_ana,cat_order_top,'Abs_TopTP')    
    fast_dotplot(adata_L1_topTP,TF_targets,TF_ana,cat_order_top,'Rel_TopTP')

adata_L1_topTP_CDay28 = adata_L1_topTP[adata_L1_topTP.obs['timePoint'] == 'Day28'].copy()
genes = TFs_dict_L1['L1_I_CD28']

## Calculate the mean expression for each gene across all cells
mean_expression = adata_L1_topTP_CDay28[:, genes].X.mean(axis=0)

# Calculate the percentage of cells expressing each gene (non-zero values)
percent_cells_expressing = (adata_L1_topTP_CDay28[:, genes].X > 0).mean(axis=0)

# Create a DataFrame to store the results
expression_df = pd.DataFrame({
    "Gene": genes,
    "MeanExpression": mean_expression.flatten(),  # Convert sparse matrix to array
    "PercentCellsExpressing": percent_cells_expressing.flatten()  # Sparse to array
})

# Filter the DataFrame for genes with mean expression > 0.5 and percent cells expressing > 20%
filtered_genes = expression_df[(expression_df['MeanExpression'] > 0) & 
                               (expression_df['PercentCellsExpressing'] > 0.02)]

# Print the filtered DataFrame
print(filtered_genes)

TF_targets = filtered_genes['Gene'].to_list()
fast_dotplot(adata_L1,TF_targets,TF_ana,cat_order_all,'Rel_AllTP_expressed')
fast_dotplot(adata_L1_topTP,TF_targets,TF_ana,cat_order_top,'Rel_TopTP_expressed')



#L2
adata_L2 = adata[np.isin(adata.obs['Annotation'],['L2', 'Ambig_L1|L2'])].copy()
adata_L2_topTP = adata_L2[np.isin(adata_L2.obs['timePoint'],cat_order_top.tolist())].copy()

for TF_ana in TFs_dict_L2.keys():
    print(TF_ana)
    
    if not TFs_dict_L2[TF_ana]:
        print(f'{TF_ana} is empty')
        continue

    TF_targets = [g for g in TFs_dict_L2[TF_ana] if g in adata_L2.var.index]
    fast_dotplot(adata_L2,TF_targets,TF_ana,cat_order_all,'Abs_AllTP')
    fast_dotplot(adata_L2,TF_targets,TF_ana,cat_order_all,'Rel_AllTP')
    
    TF_targets = [g for g in TFs_dict_L2[TF_ana] if g in adata_L2_topTP.var.index]
    fast_dotplot(adata_L2_topTP,TF_targets,TF_ana,cat_order_top,'Abs_TopTP')    
    fast_dotplot(adata_L2_topTP,TF_targets,TF_ana,cat_order_top,'Rel_TopTP')

# Subset per Luminal Types and timepoints used for TF analysis
#L1
for TF_ana in TFs_dict_L1.keys():
    print(TF_ana)
    info = TF_ana.split('_')

    TF_targets = [g for g in TFs_dict_L1[TF_ana] if g in adata_L1.var.index]

    if (info[1] == 'I') & (info[2] == 'CD28'):
        adata_inUse = adata_L1[np.isin(adata_L1.obs['timePoint'],['Intact','Day28'])].copy()
        cat_order = pd.Categorical([ 'Intact','Day28'], categories=['Intact','Day28'])
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Abs_IvsCD28')    
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Rel_IvsCD28')    


    elif (info[1] == 'CD28') & (info[2] == 'RD3'):
        adata_inUse = adata_L1[np.isin(adata_L1.obs['timePoint'],['RegenDay3','Day28'])].copy()
        cat_order = pd.Categorical(['Day28','RegenDay3'], categories=['Day28','RegenDay3'])
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Abs_CD28vsRD3')    
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Rel_CD28vsRD3')   

#L2
for TF_ana in TFs_dict_L2.keys():
    print(TF_ana)
    if not TFs_dict_L2[TF_ana]:
        print(f'{TF_ana} is empty')
        continue

    info = TF_ana.split('_')
    TF_targets = [g for g in TFs_dict_L2[TF_ana] if g in adata_L2.var.index]

    if (info[1] == 'I') & (info[2] == 'CD28'):
        adata_inUse = adata_L2[np.isin(adata_L2.obs['timePoint'],['Intact','Day28'])].copy()
        cat_order = pd.Categorical([ 'Intact','Day28'], categories=['Intact','Day28'])
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Abs_IvsCD28')    
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Rel_IvsCD28')  

    elif (info[1] == 'CD28') & (info[2] == 'RD3'):
        adata_inUse = adata_L2[np.isin(adata_L2.obs['timePoint'],['RegenDay3','Day28'])].copy()
        cat_order = pd.Categorical(['Day28','RegenDay3'], categories=['Day28','RegenDay3'])
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Abs_CD28vsRD3')    
        fast_dotplot(adata_inUse,TF_targets,TF_ana,cat_order,'Rel_CD28vsRD3')   


############################
#### TFs separetely, per TimePoints per annot
############################


adata_topTP = adata[np.isin(adata.obs['timePoint'],cat_order_top.tolist())].copy()

adata.obs['Annotation']=adata.obs['Annotation'].cat.remove_unused_categories()
# check if a gene is present
#adata.var_names[adata.var_names == 'Krt8']

for tf_dic in [TFs_dict_L1,TFs_dict_L2,{'Sanity':['Actb','Epcam','Krt8']}]:
    for TF_ana in tf_dic.keys():
        print('\n\n',TF_ana)
        
        if not tf_dic[TF_ana]:
            print(f'{TF_ana} is empty')
            continue

        TF_targets = [g for g in tf_dic[TF_ana] if g in adata_topTP.var.index]

        for gene_name in TF_targets:
            print(f'DotPlotting {TF_ana} {gene_name}')
            annot_DotPlot(adata_topTP,TF_ana,cat_order_top,gene_name)
            annot_DotPlot(adata,TF_ana,cat_order_all,gene_name)

GoI = ['Fos','Fosb','Fosl1','Fosl2','Jun','Junb']
TF_ana = 'Req'
for gene_name in GoI:
    print(f'DotPlotting{gene_name}')
    annot_DotPlot(adata_topTP,TF_ana,cat_order_top,gene_name)
    annot_DotPlot(adata,TF_ana,cat_order_all,gene_name)

TF_targets = [g for g in TopHits_250219 if g in adata_L1_topTP.var.index]
TF_ana = 'TopHits_250219'
for gene_name in TF_targets:
    print(f'DotPlotting{gene_name}')
    annot_DotPlot(adata_L1_topTP,TF_ana,cat_order_top,gene_name)
    annot_DotPlot(adata,TF_ana,cat_order_all,gene_name)

############################
#### STAT SIGNI FOR I VS CD28 
############################

TF_ana = 'L1_I_CD28'
adata_I_CD28 = adata[np.isin(adata.obs['timePoint'],['Intact','Day28'])].copy()

TF_targets = [g for g in TFs_dict_L1[TF_ana] if g in adata_I_CD28.var.index]
cat_order_I_CD28 = pd.Categorical([ 'Intact','Day28'], categories=['Intact','Day28'])

sc.tl.rank_genes_groups(adata_I_CD28, groupby='timePoint', groups=[ 'Intact','Day28'],reference='rest', method="wilcoxon",use_raw=False,pts=True)

result = adata_I_CD28.uns['rank_genes_groups'].copy()
groups = result['names'].dtype.names
df1 = pd.DataFrame({group+'_' + key:result[key][group] for group in groups for key in ['names','scores','logfoldchanges','pvals','pvals_adj']})
df2 = pd.DataFrame({group+'_' + key:result[key][group] for group in groups for key in ['pts','pts_rest']})
#to_csv
#pd.concat([df1[[group+'_names',group+'_scores',group+'_logfoldchanges',group+'_pvals',group+'_pvals_adj']].merge(df2[[group+"_pts",group+"_pts_rest"]],how="left",left_on=group+"_names",right_index=True) for group in groups],axis=1).to_csv("IvsCD28_DEG.csv")

output_file = "IvsCD28_DEG.xlsx"
with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    # Write params as a separate sheet
    params_df = pd.DataFrame(list(result['params'].items()), columns=['Parameter', 'Value'])
    params_df.to_excel(writer, sheet_name='Params', index=False)

    # Write pts data to a sheet
    df2.to_excel(writer, sheet_name='Fraction cells expr genes', index=True)

    # Write names, scores, pvals, pvals_adj, logfoldchanges to another sheet
    df1.to_excel(writer, sheet_name='Stats', index=False)

for valll in ['scores', 'logfoldchanges', 'pvals_adj']:
    sc.pl.dotplot.DEFAULT_LEGENDS_WIDTH = 4  # inches
    if valll=='logfoldchanges':
        pl_fig = sc.pl.rank_genes_groups_dotplot(adata_I_CD28,var_names=TF_targets,
            groups='Day28',
            values_to_plot=valll,
            categories_order=cat_order_I_CD28,
            cmap='bwr',
            vmin=-3,
            vmax=3,
            swap_axes=True,
            return_fig=True
            )
        pl_fig.add_totals(color='lightpink').savefig(f'./figures/{Experiment}_{TF_ana}_{valll}.pdf')
        plt.close()
    elif valll=='pvals_adj':
        pl_fig = sc.pl.rank_genes_groups_dotplot(adata_I_CD28,var_names=TF_targets,
            groups='Day28',
            values_to_plot=valll,
            categories_order=cat_order_I_CD28,
            vmin=0,
            vmax=0.1,
            swap_axes=True,
            return_fig=True
            )
        pl_fig.add_totals(color='lightpink').savefig(f'./figures/{Experiment}_{TF_ana}_{valll}.pdf')
        plt.close()
    else:
        pl_fig = sc.pl.rank_genes_groups_dotplot(adata_I_CD28,var_names=TF_targets,
            groups='Day28',
            values_to_plot=valll,
            categories_order=cat_order_I_CD28,
            swap_axes=True,
            return_fig=True
            )
        pl_fig.add_totals(color='lightpink').savefig(f'./figures/{Experiment}_{TF_ana}_{valll}.pdf')
        plt.close()




############################
#### WOUTER
############################

Experiment='Wouter2020'

adata = sc.read_h5ad(adata_dir + 'Wouter2020.h5ad').copy()
# Load the metadata TSV file
metadata = pd.read_csv("Wouter2020RawFiles/GSE146811_mmProstate10x_full_sample_final.tsv", sep="\t")

# Ensure sampleID is unique (no duplicates)
if metadata["sampleID"].duplicated().any():
    raise ValueError("Duplicate sampleID entries found in the metadata file!")

# Check if all sampleIDs exist in adata.obs
adata_cells = set(adata.obs.index)
metadata_cells = set(metadata["sampleID"])
missing_in_adata = metadata_cells - adata_cells
missing_in_metadata = adata_cells - metadata_cells
if missing_in_adata:
    raise ValueError(f"The following sampleIDs from metadata are missing in adata: {missing_in_adata}")
if missing_in_metadata:
    raise ValueError(f"The following sampleIDs from adata are missing in metadata: {missing_in_metadata}")

# Set index for merging
metadata.set_index("sampleID", inplace=True)
# Ensure 1:1 matching by reordering metadata to match adata.obs
metadata = metadata.loc[adata.obs.index]
# Add metadata columns to adata.obs
adata.obs = adata.obs.join(metadata)
adata.obs['Annotation'] = adata.obs['predType'].copy()
# Convert columns to categorical
adata.obs["timePoint"] = adata.obs["timePoint"].astype("category")
adata.obs["Annotation"] = adata.obs["Annotation"].astype("category")
adata.X = adata.X.toarray()

adata = adata[adata.obs['Annotation'].str.startswith(('Epi_Luminal','Epi_Basal'))].copy()
adata.obs['Annotation'].cat.remove_unused_categories()

cat_order_all = pd.Categorical(['intact1','intact2','CastDay1', 'CastDay7', 'CastDay14', 'CastDay28', 'RegenDay1', 'RegenDay2', 'RegenDay3','RegenDay7', 'RegenDay14', 'RegenDay28'],categories=['intact1','intact2','CastDay1', 'CastDay7', 'CastDay14', 'CastDay28', 'RegenDay1', 'RegenDay2', 'RegenDay3','RegenDay7', 'RegenDay14', 'RegenDay28'])
adata = adata[np.isin(adata.obs['timePoint'],cat_order_all)].copy()
adata.obs['timePoint'].cat.remove_unused_categories()

cat_order_top = pd.Categorical(['intact1','intact2','CastDay28'],categories=['intact1','intact2','CastDay28'])
adata_TP = adata[np.isin(adata.obs['timePoint'],cat_order_top)].copy()
adata_TP.obs['timePoint'].cat.remove_unused_categories()

TF_targets = [g for g in TopHits_250219 if g in adata.var.index]
TF_ana = 'TopHits_250219'
for gene_name in TF_targets:
    print(f'DotPlotting {gene_name}')
    annot_DotPlot(adata_TP,TF_ana,cat_order_top,gene_name,Experiment)
    annot_DotPlot(adata,TF_ana,cat_order_all,gene_name,Experiment)

adata_TP_L1 = adata_TP[np.isin(adata_TP.obs['Annotation'],'Epi_Luminal_1')].copy()
adata_TP_L1.obs['Annotation'].cat.remove_unused_categories()

adata_L1 = adata[np.isin(adata.obs['Annotation'],'Epi_Luminal_1')].copy()
adata_L1.obs['Annotation'].cat.remove_unused_categories()

TF_targets = [g for g in TopHits_250219 if g in adata.var.index]
TF_ana = 'TopHits_250219'
fast_dotplot(adata,TF_targets,TF_ana,cat_order_all,'Abs_AllTP',Experiment) 
fast_dotplot(adata,TF_targets,TF_ana,cat_order_all,'Rel_AllTP',Experiment) 

fast_dotplot(adata_L1,TF_targets,TF_ana,cat_order_all,'Abs_AllTP_L1',Experiment) 
fast_dotplot(adata_L1,TF_targets,TF_ana,cat_order_all,'Rel_AllTP_L1',Experiment) 

fast_dotplot(adata_TP,TF_targets,TF_ana,cat_order_top,'Abs_AllTPOnly',Experiment) 
fast_dotplot(adata_TP,TF_targets,TF_ana,cat_order_top,'Rel_AllTPOnly',Experiment) 

fast_dotplot(adata_TP_L1,TF_targets,TF_ana,cat_order_top,'Abs_AllTPOnlyL1',Experiment) 
fast_dotplot(adata_TP_L1,TF_targets,TF_ana,cat_order_top,'Rel_AllTPOnlyL1',Experiment) 



TF_ana = 'L1_I_CD28'
TF_targets = [g for g in TFs_dict_L1[TF_ana] if g in adata.var.index]

fast_dotplot(adata,TF_targets,TF_ana,cat_order_all,'Abs_AllTP',Experiment) 
fast_dotplot(adata,TF_targets,TF_ana,cat_order_all,'Rel_AllTP',Experiment) 

cat_order_all = pd.Categorical(['intact1','intact2','CastDay28','RegenDay3'],categories=['intact1','intact2','CastDay28','RegenDay3'])
adata_TP = adata[np.isin(adata.obs['timePoint'],cat_order_all)]
adata_TP.obs['timePoint'].cat.remove_unused_categories()
fast_dotplot(adata_TP,TF_targets,TF_ana,cat_order_all,'Abs_AllTPOnly',Experiment) 
fast_dotplot(adata_TP,TF_targets,TF_ana,cat_order_all,'Rel_AllTPOnly',Experiment) 














































































def dotplot(TF_targets,TFs_dict=TFs_dict,adata=adata,TP_cat_order=cat_order):

        print(f"Dotplotting {TF_name}")

        if TF_name in adata.var.index:
            sc.pl.dotplot(adata, var_names=TF_name, groupby='timePoint', use_raw=False,
                standard_scale='var',cmap='bwr',swap_axes=True,categories_order=cat_order,# dot_max=0.5,
                title=f'{Experiment} \n{TF_name} Solo', save=f'{Experiment}_Solo_{TF_name}_{timepointUsed}_dotplot', show=False)
        else:
            print(f"{TF_name} itself not expressed")

        if df[TF_name].dropna().empty:
            print(TF_name+' is empty')
            continue

        df_targets = df[[TF_name,TF_name+'_mor']]
        df_targets = df_targets.dropna()
        TF_targets = df_targets[TF_name]

        tar_present_in_adata = [gene for gene in TF_targets if gene in adata.var.index]
        tar_absent_in_adata = [gene for gene in TF_targets if gene not in adata.var.index]
        print(f"{TF_name} regulated genes were absent:",tar_absent_in_adata)
        print(f"{TF_name} regulated genes are present:",tar_present_in_adata,'\n')

        df_targets = df_targets[df_targets[TF_name].isin(tar_present_in_adata)]
        TF_targets = df_targets[TF_name]
        # Sorting the DataFrame based on the second column
        df_sorted = df_targets.sort_values(by=TF_name+'_mor').reset_index(drop=True)

        # Initializing lists to store positions and corresponding values
        positions = []
        values = []

        # Identifying consecutive same values
        start_idx = 0

        for i in range(1, len(df_sorted)):
            if df_sorted[f'{TF_name}_mor'][i] != df_sorted[f'{TF_name}_mor'][start_idx]:
                if i - start_idx > 1:  # There are consecutive values
                    positions.append((start_idx, i - 1))
                    values.append(round(df_sorted[f'{TF_name}_mor'][start_idx], 1))
                start_idx = i

        # Checking the last group of consecutive values
        if len(df_sorted) - start_idx > 1:
            positions.append((start_idx, len(df_sorted) - 1))
            values.append(round(df_sorted[f'{TF_name}_mor'][start_idx], 1))

        # Handling single entries (not consecutive but unique values)
        for i in range(len(df_sorted)):
            if i == 0 or df_sorted[f'{TF_name}_mor'][i] != df_sorted[f'{TF_name}_mor'][i - 1]:
                if i == len(df_sorted) - 1 or df_sorted[f'{TF_name}_mor'][i] != df_sorted[f'{TF_name}_mor'][i + 1]:
                    positions.append((i, i))
                    values.append(round(df_sorted[f'{TF_name}_mor'][i], 1))

        # Sorting positions and values to ensure correct order after adding single entries
        sorted_positions_values = sorted(zip(positions, values))
        positions, values = zip(*sorted_positions_values)

        positions = list(positions)
        values = [str(va) for va in values]

        if len(TF_targets) >120:
            print(f"Too many genes regulated to plot {TF_name} regulated genes dotplot")
            continue

        sc.pl.dotplot(adata, var_names=TF_targets, groupby='timePoint', use_raw=False,
            standard_scale='var',cmap='bwr',swap_axes=True,categories_order=cat_order,# dot_max=0.5,
            var_group_positions=positions,var_group_labels=values,
            title=f'{Experiment} {db_name}\n{TF_name}', save=f'{Experiment}_{db_name}_{TF_name}_{timepointUsed}_dotplot', show=False)



#Change if other celltype
#adata = adata[adata.obs['TypeDiffAmbig']=='L1 24923cells']
adata = adata[np.isin(adata.obs['TypeDiffAmbig'],['Ambig_L1|L2 845cells', 'L2 559cells'])]

#Change if other timepoints
#adata = adata[np.isin(adata.obs["timePoint"],['Day28', 'Intact'])]
adata = adata[np.isin(adata.obs["timePoint"],['Day28','RegenDay3'])]
#cat_order = pd.Categorical(['Intact','Day28'], categories=['Intact','Day28'])
cat_order = pd.Categorical(['Day28','RegenDay3'], categories=['Day28','RegenDay3'])
#cat_order = pd.Categorical([ 'Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'], categories=[ 'Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'])


timepointUsed = "IntactvsDay28"#"Day28vsRegenDay3"# 'AllTPs', 
#adata import
if 'SING' in Experiment:
    adata = adata[~np.isin(adata.obs["batch"],['WK-1501_BL6_INTACT_AP_Test3_SORT', 'WK-1384_BL6_Intact_AP_2_SLT'])]
    
    #Change if other celltype
    #adata = adata[adata.obs['TypeDiffAmbig']=='L1 24923cells']
    adata = adata[np.isin(adata.obs['TypeDiffAmbig'],['Ambig_L1|L2 845cells', 'L2 559cells'])]
    
    #Change if other timepoints
    #adata = adata[np.isin(adata.obs["timePoint"],['Day28', 'Intact'])]
    adata = adata[np.isin(adata.obs["timePoint"],['Day28','RegenDay3'])]
    #cat_order = pd.Categorical(['Intact','Day28'], categories=['Intact','Day28'])
    cat_order = pd.Categorical(['Day28','RegenDay3'], categories=['Day28','RegenDay3'])
    #cat_order = pd.Categorical([ 'Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'], categories=[ 'Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'])
if 'ENZ' in Experiment:
    adata = adata[adata.obs['tissueProv'].str.startswith('AP')]
    adata = adata[adata.obs['TypeDiffAmbig']=='L1 7108cells']
    cat_order = pd.Categorical(['Intact','RegenDay3'], categories=['Intact','RegenDay3'])


FGF7_targets = ["Vegfa","Sftpd","Wnt4","Ccdc3"]
FGF2_targets = ["Hes1","Cdkn1a","Fosl1","Tgfb1","Stc1","Ier3","Ltbp1","Cdcp1","Pmaip1","Tmem158"]
TGFb_targets = ["Serpine1","Edn1","Angptl4","Smad7","Id3","Gadd45b","Klf10","Sphk1","Tgfbr1","Junb"]
IL6_targets = ["Socs3","Stat3","Junb","Ccnd1","Hif1a","Icam1","Bcl6","Ass1"]

all_targets = [(FGF7_targets,"FGF7_Targets"),
(FGF2_targets,"FGF-2_Targets"),
(TGFb_targets,"TGF-b_Targets"),
(IL6_targets,"IL6_Targets")]

for tar,nam in all_targets:
    print(f"Dotplotting {nam}")

    tar_present_in_adata = [gene for gene in tar if gene in adata.var.index]
    tar_absent_in_adata = [gene for gene in tar if gene not in adata.var.index]
    print("These genes were absent:",tar_absent_in_adata)
    print("These genes are present:",tar_present_in_adata,'\n')

    sc.pl.dotplot(adata, var_names=tar_present_in_adata, groupby='timePoint', use_raw=False,
        standard_scale='var',cmap='bwr',swap_axes=True,categories_order=cat_order,# dot_max=0.5,
        title=f'{Experiment}\n{nam}', save=f'{Experiment}_{nam}_{timepointUsed}_dotplot', show=False)

def dotplot(pathtsv,db_name,adata=adata,cat_order=cat_order,timepointUsed=timepointUsed):
    df = pd.read_csv(pathtsv, sep='\t',header=0, index_col=False)

    for TF_name in df.columns:
        if ('mor' in TF_name):
            continue
        print(f"Dotplotting {TF_name}")

        if TF_name in adata.var.index:
            sc.pl.dotplot(adata, var_names=TF_name, groupby='timePoint', use_raw=False,
                standard_scale='var',cmap='bwr',swap_axes=True,categories_order=cat_order,# dot_max=0.5,
                title=f'{Experiment} \n{TF_name} Solo', save=f'{Experiment}_Solo_{TF_name}_{timepointUsed}_dotplot', show=False)
        else:
            print(f"{TF_name} itself not expressed")

        if df[TF_name].dropna().empty:
            print(TF_name+' is empty')
            continue

        df_targets = df[[TF_name,TF_name+'_mor']]
        df_targets = df_targets.dropna()
        TF_targets = df_targets[TF_name]

        tar_present_in_adata = [gene for gene in TF_targets if gene in adata.var.index]
        tar_absent_in_adata = [gene for gene in TF_targets if gene not in adata.var.index]
        print(f"{TF_name} regulated genes were absent:",tar_absent_in_adata)
        print(f"{TF_name} regulated genes are present:",tar_present_in_adata,'\n')

        df_targets = df_targets[df_targets[TF_name].isin(tar_present_in_adata)]
        TF_targets = df_targets[TF_name]
        # Sorting the DataFrame based on the second column
        df_sorted = df_targets.sort_values(by=TF_name+'_mor').reset_index(drop=True)

        # Initializing lists to store positions and corresponding values
        positions = []
        values = []

        # Identifying consecutive same values
        start_idx = 0

        for i in range(1, len(df_sorted)):
            if df_sorted[f'{TF_name}_mor'][i] != df_sorted[f'{TF_name}_mor'][start_idx]:
                if i - start_idx > 1:  # There are consecutive values
                    positions.append((start_idx, i - 1))
                    values.append(round(df_sorted[f'{TF_name}_mor'][start_idx], 1))
                start_idx = i

        # Checking the last group of consecutive values
        if len(df_sorted) - start_idx > 1:
            positions.append((start_idx, len(df_sorted) - 1))
            values.append(round(df_sorted[f'{TF_name}_mor'][start_idx], 1))

        # Handling single entries (not consecutive but unique values)
        for i in range(len(df_sorted)):
            if i == 0 or df_sorted[f'{TF_name}_mor'][i] != df_sorted[f'{TF_name}_mor'][i - 1]:
                if i == len(df_sorted) - 1 or df_sorted[f'{TF_name}_mor'][i] != df_sorted[f'{TF_name}_mor'][i + 1]:
                    positions.append((i, i))
                    values.append(round(df_sorted[f'{TF_name}_mor'][i], 1))

        # Sorting positions and values to ensure correct order after adding single entries
        sorted_positions_values = sorted(zip(positions, values))
        positions, values = zip(*sorted_positions_values)

        positions = list(positions)
        values = [str(va) for va in values]

        if len(TF_targets) >120:
            print(f"Too many genes regulated to plot {TF_name} regulated genes dotplot")
            continue

        sc.pl.dotplot(adata, var_names=TF_targets, groupby='timePoint', use_raw=False,
            standard_scale='var',cmap='bwr',swap_axes=True,categories_order=cat_order,# dot_max=0.5,
            var_group_positions=positions,var_group_labels=values,
            title=f'{Experiment} {db_name}\n{TF_name}', save=f'{Experiment}_{db_name}_{TF_name}_{timepointUsed}_dotplot', show=False)


if 'SING' in Experiment:
    dotplot("HOCOMOCOv12_L1_IntactvDay28.tsv","HOCOMOCOv12")
    dotplot("JASPAR2022CORE_L1_IntactvDay28.tsv","JASPAR2022CORE")
#dotplot("JASPAR2022CORE_L2_AmbigL1L2_IntactvDay28.tsv","JASPAR2022CORE")

dotplot("HOCOMOCOv12_L2A_IntactvDay28.tsv","HOCOMOCOv12")

def dotplot_cols(pathtsv,adata=adata,timepointUsed=timepointUsed):
    df = pd.read_csv(pathtsv, sep='\t',header=0, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('_mor')]
    TF_targets = df.columns.tolist()

    tar_present_in_adata = [gene for gene in TF_targets if gene in adata.var.index]
    tar_absent_in_adata = [gene for gene in TF_targets if gene not in adata.var.index]
    print(f"{pathtsv} regulated genes were absent:",tar_absent_in_adata)
    print(f"{pathtsv} regulated genes are present:",tar_present_in_adata,'\n')

    sc.pl.dotplot(adata, var_names=tar_present_in_adata, groupby='timePoint', use_raw=False,
        standard_scale='var',cmap='bwr',swap_axes=True,categories_order=cat_order,# dot_max=0.5,
        title=f'{Experiment} \n{pathtsv.rstrip(".tsv")}', save=f'{Experiment}_TFExpression_{timepointUsed}_dotplot', show=False)

dotplot_cols("HOCOMOCOv12_L1_IntactvDay28.tsv")
dotplot_cols("HOCOMOCOv12_L2A_IntactvDay28.tsv")

dotplot_cols("HOCOMOCOv12_L1_Day28vRegenDay3.tsv")
dotplot_cols("HOCOMOCOv12_L2A_Day28vRegenDay3.tsv")


