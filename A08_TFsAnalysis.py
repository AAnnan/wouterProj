##Downstream Analysis

# Purpose:
#   define open chromatin and TF’s binding motifs enriched
#   enriched in one vs all 

#   for example in the singular: enriched in INTACT VS. D28/RD3
#   or pairwise comparision (intact vs. castrationD28 and intact vs. RD3 etc)

import snapatac2 as snap
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import polars as pl
import os
import warnings
from collections import Counter
warnings.simplefilter(action='ignore', category=FutureWarning)
sc.settings.verbosity = 0
sc.settings.set_figure_params(
	figsize=(6, 6),
    dpi_save=300,
    fontsize=12,
    facecolor="white",
    frameon=True
    ,
)
print(f'Anndata: {ad.__version__}\nSnapatac2: {snap.__version__}\nScanpy: {sc.__version__}')

# To change according to samples
# SING for singulator, ENZ for enzymatic digestion
Experiment='Wouter21_SING'
temp_macs3 = '/mnt/etemp/ahrmad/wouter/Wouter21_SING/annot/temp'
TFmotifs = '/mnt/etemp/ahrmad/wouter/motif_databases/CIS-BP_2.00/Mus_musculus.meme'
TFmotifs = '/mnt/etemp/ahrmad/wouter/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme'
TFmotifs = '/mnt/etemp/ahrmad/wouter/motif_databases/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme'
TFmotifs = '/mnt/etemp/ahrmad/wouter/motif_databases/H12CORE_meme_format.meme'

atac_ext = '_Post.h5ad'
rna_h5ad = "/mnt/etemp/ahrmad/wouter/RNA_CB/annot_other_norm/Wouter21_SING_CB_annot_All_newSigs.h5ad"

leiR = 1
leiRes = f'leiden_res{leiR}'

adata_ATAC = snap.read(Experiment + atac_ext, backed=None)
adata_AllAnnot = sc.read_h5ad(rna_h5ad)

# Add RNA All Annot to adata_ATAC
adata_ATAC = adata_ATAC[np.isin(adata_ATAC.obs.index,adata_AllAnnot.obs.index)]
adata_ATAC.obs = adata_ATAC.obs.reindex(index=sorted(adata_ATAC.obs.index))
adata_AllAnnot.obs = adata_AllAnnot.obs.reindex(index=sorted(adata_AllAnnot.obs.index))

assert (adata_ATAC.obs.index == adata_AllAnnot.obs.index).all()
adata_ATAC.obs['Annotation'] = adata_AllAnnot.obs['Annotation']
if np.sum(adata_ATAC.obs['Annotation'] == 'Other')==0:
    adata_ATAC.obs['Annotation'] = adata_ATAC.obs['Annotation'].cat.remove_categories('Other')

desired_order = ['Intact', 'Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3']
adata_ATAC.obs['timePoint'] = adata_ATAC.obs['timePoint'].astype(pd.CategoricalDtype(categories=desired_order, ordered=True))

# Plotting
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['Annotation'],save=Experiment+'_TF_Annot',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['timePoint'],save=Experiment+'_TF_timePoints',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['batch'],save=Experiment+'_TF_batch',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['batch'],save=Experiment+'_TF_batch',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=["n_fragment", "frac_dup", "frac_mito"],save=Experiment+'_Counts_QC',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=["FRiP", "tsse"],save=Experiment+'_PctCounts_QC',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=leiRes,legend_loc="on data",save=f"{Experiment}_cluster{leiRes}",title=Experiment,show=False)

#for t in ['batch','timePoint','Annotation']:
#    print('umap',t)
#    snap.pl.umap(adata_ATAC,marker_size=75000/adata_ATAC.n_obs,marker_opacity=0.6,color=t,out_file=f'figures/umap{Experiment}_TF_QC_SNAP_{t}.pdf',show=False)

wanted_tp_comps = [['Intact','Day28'],['Day28','RegenDay3']] #,('Day28','RegenDay1'),('Day28','RegenDay2')

#Subset by Annotation
luminal_celltypes = [['L1'], ['L2', 'Ambig_L1|L2']]#,['Basal']]

for l_types in luminal_celltypes:
    if l_types == ['L1']:
        log_fc_min = 1
    elif l_types == ['L2', 'Ambig_L1|L2']:
        log_fc_min = 0.3

    adata_ATAC_subset = adata_ATAC[np.isin(adata_ATAC.obs['Annotation'],l_types)].copy()
    #adata_ATAC_subset = adata_ATAC_subset[adata_ATAC_subset.obs['batch']!='WK-1585_INTACT_AP_BL6_Citrate']
    print(Counter(adata_ATAC_subset.obs['batch']))

    # Finding marker regions
    # aggregates signal across cells and utilizes z-scores to identify specifically enriched peaks.
    # Does not consider the variations across cells.
    for TPcomp in wanted_tp_comps:
        TPcomp_list = TPcomp
        print(TPcomp)
        adata_ATAC_subset_comp = adata_ATAC_subset[np.isin(adata_ATAC_subset.obs['timePoint'],TPcomp)].copy()
        print(Counter(adata_ATAC_subset_comp.obs['batch']))
        #Annotation
        marker_peaks = snap.tl.marker_regions(adata_ATAC_subset_comp, groupby='timePoint', pvalue=0.16)

        # Resulting heatmap
        snap.pl.regions(adata_ATAC_subset_comp, groupby='timePoint', width=1200, height=800, peaks=marker_peaks, out_file=f'figures/heatmap{Experiment}_TF_{l_types}_regionsTP_{"_".join(TPcomp_list)}.png',show=False)

        TFmotifs = '/mnt/etemp/ahrmad/wouter/motif_databases/H12CORE_meme_format.meme'

        motifs = snap.tl.motif_enrichment(
            motifs=snap.read_motifs(TFmotifs),
            regions=marker_peaks,
            genome_fasta=snap.genome.mm10,
        )
        
        snap.pl.motif_enrichment(motifs, max_fdr=0.001, min_log_fc=log_fc_min, height=1500, out_file=f'figures/heatmap{Experiment}_TF_motifsHOCOMOCOv12_{l_types}_TP_{"_".join(TPcomp_list)}.png',show=False)
        #snap.pl.motif_enrichment(motifs, max_fdr=0.001, min_log_fc=1, height=1500, out_file=f'figures/heatmap{Experiment}_TF_motifsHOCOMOCOv12_{l_types}_TP_{"_".join(TPcomp_list)}.png',show=False)
        #snap.pl.motif_enrichment(motifs, max_fdr=0.001, min_log_fc=1, height=1500, out_file=f'figures/heatmap{Experiment}_TF_motifsHOCOMOCOv12_{l_types}_TP_{"_".join(TPcomp_list)}.html',show=False, interactive=True)
        
        TFmotifs = '/mnt/etemp/ahrmad/wouter/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme'

        motifs = snap.tl.motif_enrichment(
            motifs=snap.read_motifs(TFmotifs),
            regions=marker_peaks,
            genome_fasta=snap.genome.mm10,
        )

        snap.pl.motif_enrichment(motifs, max_fdr=0.001, min_log_fc=log_fc_min, height=1500, out_file=f'figures/heatmap{Experiment}_TF_motifsJASPAR2022CORE_{l_types}_TP_{"_".join(TPcomp_list)}.png',show=False)
        #snap.pl.motif_enrichment(motifs, max_fdr=0.001, min_log_fc=1, height=1500, out_file=f'figures/heatmap{Experiment}_TF_motifsJASPAR2022CORE_{l_types}_TP_{"_".join(TPcomp_list)}.png',show=False)
        #snap.pl.motif_enrichment(motifs, max_fdr=0.001, min_log_fc=1, height=1500, out_file=f'figures/heatmap{Experiment}_TF_motifsJASPAR2022CORE_{l_types}_TP_{"_".join(TPcomp_list)}.html',show=False, interactive=True)

        # AHR.H12CORE.0.P.B as an example.
        #Here AHR denotes the UniProt AC prefix (most of the time identical between human and mouse orthologs, e.g. AHR_HUMAN and AHR_MOUSE).
        #H12CORE denotes the subcollection, and can also be H12RSNP/H12INVIVO/H12INVITRO in downloadable motifs sets.
        #0 is the subtype number, where 0 denotes the most common motifs scoring the best across all benchmarking datasets.
        #P is the type of the experiment that yielded motifs that were assigned to the subtype during expert curation. Can be P (ChIP-Seq), S (HT-SELEX), or M (Methyl-HT-SELEX), or any combination of those three for motifs found in several types of experiments.
        #B is the motif quality on the ABCD scale, see below. 


    #Subset by Annotation

        #Peak calling at the timepoint level
        snap.tl.macs3(adata_ATAC_subset_comp, groupby='timePoint',tempdir=temp_macs3)
        peaks = snap.tl.merge_peaks(adata_ATAC_subset_comp.uns['macs3'], snap.genome.mm10)
        peaks.head()
        peak_mat = snap.pp.make_peak_matrix(adata_ATAC_subset_comp, use_rep=peaks['Peaks'])
        peak_mat


        # Regression-based differential test

        def diff_test(g1,g2,group_type='timePoint',peaks=peaks,peak_mat=peak_mat,data=adata_ATAC_subset_comp):
            g1_cells = data.obs[group_type] == g1
            g2_cells = data.obs[group_type] == g2
            peaks_selected = np.logical_or(peaks[g1].to_numpy(), peaks[g2].to_numpy())

            # Perform differential test using tl.diff_test.
            diff_peaks = snap.tl.diff_test(
                peak_mat,
                cell_group1=g1_cells,
                cell_group2=g2_cells,
                features=peaks_selected,
            )

            diff_peaks = diff_peaks.filter(pl.col('adjusted p-value') < 0.01)
            #diff_peaks.head()

            snap.pl.regions(
                peak_mat,
                groupby = group_type,
                peaks = {
                    g1: diff_peaks.filter(pl.col("log2(fold_change)") > 0)['feature name'].to_numpy(),
                    g2: diff_peaks.filter(pl.col("log2(fold_change)") < 0)['feature name'].to_numpy(),
                },
                out_file=f'figures/heatmap{Experiment}_TF_Diff_{g1}vs{g2}_{"_".join(TPcomp_list)}.png',show=False
            )
            return 0

        diff_test(TPcomp[0],TPcomp[1],group_type='timePoint')

def diff_test_vsALL(g1,group_type='timePoint',peaks=peaks,peak_mat=peak_mat,data=adata_ATAC_subset):

    g1_cells = data.obs[group_type] == g1

    barcodes = np.array(data.obs_names)
    background = []
    for i in np.unique(data.obs[group_type]):
        if i != g1:
            cells = np.random.choice(barcodes[data.obs[group_type] == i], size=50, replace=False)
            background.append(cells)
    background = np.concatenate(background)

    diff_peaks = snap.tl.diff_test(
        peak_mat,
        cell_group1=g1_cells,
        cell_group2=background,
        features=peaks[g1].to_numpy(),
        direction="positive",
    )
    
    diff_peaks = diff_peaks.filter(pl.col('adjusted p-value') < 0.01)
    diff_peaks.head()

    snap.pl.regions(
        peak_mat,
        groupby=group_type,
        peaks={ g1: diff_peaks['feature name'].to_numpy() },
        out_file=f'figures/heatmap{Experiment}_TF_Diff_{g1}vsALL.png',show=False
    )
    return 0

diff_test_vsALL("Intact",group_type='timePoint')
diff_test_vsALL("Day28",group_type='timePoint')
diff_test_vsALL("RegenDay1",group_type='timePoint')
diff_test_vsALL("RegenDay2",group_type='timePoint')
diff_test_vsALL("RegenDay3",group_type='timePoint')


