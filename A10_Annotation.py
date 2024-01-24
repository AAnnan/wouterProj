##Annotate Cells

# Purpose:
#   Cell annotation / Cell typing

import scanpy as sc
import matplotlib.pyplot as plt

# To change according to samples
# SING for singulator, ENZ for enzymatic digestion
Experiments='Wouter21_ENZ'

# Input Files
DirATAC = '/mnt/etemp/ahrmad/wouter/batch_ATAC'
Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
qc_ext = '_filt.h5ad'
refDir = '/mnt/etemp/ahrmad/wouter/refs'
resDir = '/mnt/etemp/ahrmad/wouter/'

snap_ext = '_GeneMat.h5ad'

adata = sc.read_h5ad(Experiment + snap_ext)

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






