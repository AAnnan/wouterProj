##Annotate Cells

# Purpose:
#   Cell annotation / Cell typing
import snapatac2 as snap
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import os
print('Scanpy: ',sc.__version__,'Snap: ',snap.__version__)
sc.set_figure_params(figsize=(5, 5))

refDir = '/mnt/etemp/ahrmad/wouter/refs'
### TO CHECK WHEN CHANGING SAMPLES ###
Dir = '/mnt/etemp/ahrmad/wouter/Wouter21_ENZ'
Experiment='Wouter21_ENZ'
### TO CHECK WHEN CHANGING SAMPLES ###
if not (os.path.exists(refDir)) & (os.path.exists(Dir)):
    raise ValueError('Check folder paths')
gene_set_table = f'{refDir}/table_s8_summary.txt'

post_ext = '_GeneMat.h5ad'
post_ext = '_Post.h5ad'
annot_ext = '_Annotated.h5ad'

#CREATION OF THE GENE MAT
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


### Annotation
adata = sc.read_h5ad(Experiment + post_ext)

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
    atacdf = pd.read_csv(f'{refDir}/csv/{sample}_doublet_scores_ATAC_AMU.csv')
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






