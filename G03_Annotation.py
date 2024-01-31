import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import os
sc.set_figure_params(figsize=(5, 5))


refDir = '/mnt/etemp/ahrmad/wouter/refs'
### TO CHECK WHEN CHANGING SAMPLES ###
DirRNA = '/mnt/etemp/ahrmad/wouter/ENZ_RNA'
Experiment='Wouter21_ENZ'
### TO CHECK WHEN CHANGING SAMPLES ###
if not (os.path.exists(refDir)) & (os.path.exists(DirRNA)):
    raise ValueError('Check folder paths')
gene_set_table = f'{refDir}/table_s8_summary.txt'


post_ext = '_post.h5ad'
post_ext = '_Norm_FeatSel_DimRed.h5ad'
### Annotation
adata = sc.read_h5ad(Experiment + post_ext)
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

pdfjam --nup 4x2 --landscape --a4paper figures/umap*Imm*.png --outfile SING_Imm.pdf
pdfjam --nup 4x2 --landscape --a4paper figures/umap*Str*.png --outfile SING_Str.pdf
pdfjam --nup 4x2 --landscape --a4paper figures/umap*Epi*.png --outfile SING_Epi.pdf

pdfjam --nup 4x2 --landscape --a4paper figures/umap*Imm*.png --outfile ENZ_Imm.pdf
pdfjam --nup 4x2 --landscape --a4paper figures/umap*Str*.png --outfile ENZ_Str.pdf
pdfjam --nup 4x2 --landscape --a4paper figures/umap*Epi*.png --outfile ENZ_Epi.pdf









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


