#intact: BL6 ctrl
# CD28: CD28
#mm env create -n signac
#mm install r-signac r-seurat bioconductor-bsgenome.mmusculus.ucsc.mm10 bioconductor-ensdb.mmusculus.v79 r-azimuth bioconductor-biovizbase bioconductor-scdblfinder

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(ggvenn)
library(patchwork)


# load the RNA and ATAC data
countsCD28 <- Read10X_h5("/mnt/ndata/daniele/wouter/Processed/CellRangerArc/WK-1585_Castrate_Day28_AP_BL6/outs/filtered_feature_bc_matrix.h5")
fragpathCD28 <- "/mnt/ndata/daniele/wouter/Processed/CellRangerArc/WK-1585_Castrate_Day28_AP_BL6/outs/atac_fragments.tsv.gz"

countsI <- Read10X_h5("/mnt/ndata/daniele/wouter/Processed/CellRangerArc/WK-1585_INTACT_AP_BL6_Contrl/outs/filtered_feature_bc_matrix.h5")
fragpathI <- "/mnt/ndata/daniele/wouter/Processed/CellRangerArc/WK-1585_INTACT_AP_BL6_Contrl/outs/atac_fragments.tsv.gz"

atac_countsCD28 <- countsCD28$Peaks
atac_countsI <- countsI$Peaks

# get gene annotations for hg38
grange.counts = StringToGRanges(rownames(atac_countsCD28), sep = c(":", "-"))
grange.use = seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_countsCD28 = atac_countsCD28[as.vector(grange.use), ]

# get gene annotations for hg38
grange.counts = StringToGRanges(rownames(atac_countsI), sep = c(":", "-"))
grange.use = seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_countsI = atac_countsI[as.vector(grange.use), ]

annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) = 'UCSC'
genome(annotation) = "mm10"

# create a Seurat object containing the RNA adata
CD28 <- CreateSeuratObject(
  counts = countsCD28$`Gene Expression`,
  assay = "RNA"
)

I <- CreateSeuratObject(
  counts = countsI$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
CD28[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_countsCD28,
  sep = c(":", "-"),
  fragments = fragpathCD28,
  annotation = annotation
)
I[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_countsI,
  sep = c(":", "-"),
  fragments = fragpathI,
  annotation = annotation
)

DefaultAssay(CD28) <- "ATAC"
DefaultAssay(I) <- "ATAC"

#Calculate the strength of the nucleosome signal per cell. 
#Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
CD28 <- NucleosomeSignal(CD28)
#Compute the transcription start site (TSS) enrichment score for each cell, 
#as defined by ENCODE: https://www.encodeproject.org/data-standards/terms/.
CD28 <- TSSEnrichment(CD28,region_extension=2000)
dst_scatCD28 <- DensityScatter(CD28, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(filename = "ATAC0_dst_scatCD28.pdf", plot = dst_scatCD28, width = 12, height = 6)

vln_plot <- VlnPlot(
  object = CD28,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0.2
)
ggsave(filename = "ATAC0_vln_plotCD28.pdf", plot = vln_plot, width = 12, height = 6)

I <- NucleosomeSignal(I)
I <- TSSEnrichment(I,region_extension=2000)
dst_scatI <- DensityScatter(I, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(filename = "ATAC0_dst_scatI.pdf", plot = dst_scatI, width = 12, height = 6)

vln_plot <- VlnPlot(
  object = I,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0.2
)
ggsave(filename = "ATAC0_vln_plotI.pdf", plot = vln_plot, width = 12, height = 6)


DefaultAssay(CD28) <- "RNA"
DefaultAssay(I) <- "RNA"
CD28[["percent.mt"]] <- PercentageFeatureSet(CD28, pattern = "^mt-")
CD28[["percent.ribo"]] <- PercentageFeatureSet(CD28, pattern = "^Rps|^Rpl")
CD28[["percent.hb"]] <- PercentageFeatureSet(CD28, pattern = "^Hb")
I[["percent.mt"]] <- PercentageFeatureSet(I, pattern = "^mt-")
I[["percent.ribo"]] <- PercentageFeatureSet(I, pattern = "^Rps|^Rpl")
I[["percent.hb"]] <- PercentageFeatureSet(I, pattern = "^Hb")

# Visualize RNA QC metrics as a violin plot
vln_plot <- VlnPlot(CD28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 5)
ggsave(filename = "RNA0_vln_plotCD28.pdf", plot = vln_plot, width = 12, height = 6)

vln_plot <- VlnPlot(I, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 5)
ggsave(filename = "RNA0_vln_plotI.pdf", plot = vln_plot, width = 12, height = 6)

scat_plot <- FeatureScatter(I, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = "RNA0_scat_plotI.pdf", plot = scat_plot, width = 12, height = 6)

scat_plot <- FeatureScatter(CD28, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = "RNA0_scat_plotCD28.pdf", plot = scat_plot, width = 12, height = 6)

DefaultAssay(CD28) <- "ATAC"
DefaultAssay(I) <- "ATAC"

# filter out low quality cells
CD28 <- subset(
  x = CD28,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 50000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 100 &
    nucleosome_signal < 2 &
    TSS.enrichment > 5 &
    percent.mt < 3 &
    percent.ribo < 5
)
CD28

I <- subset(
  x = I,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 50000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 100 &
    nucleosome_signal < 2 &
    TSS.enrichment > 5 &
    percent.mt < 3 &
    percent.ribo < 5
)
I

#Doublet removal 
bp <- MulticoreParam(20, RNGseed=42)
DefaultAssay(CD28) <- "RNA"
DefaultAssay(I) <- "RNA"
q <- scDblFinder(CD28$RNA$counts, BPPARAM=bp)
CD28[['scDblFinder.score.RNA']] <- q$scDblFinder.score
CD28[['scDblFinder.class.RNA']] <- q$scDblFinder.class

q <- scDblFinder(I$RNA$counts, BPPARAM=bp)
I[['scDblFinder.score.RNA']] <- q$scDblFinder.score
I[['scDblFinder.class.RNA']] <- q$scDblFinder.class

DefaultAssay(CD28) <- "ATAC"
DefaultAssay(I) <- "ATAC"
q <- scDblFinder(CD28$ATAC$counts, BPPARAM=bp)
CD28[['scDblFinder.score.ATAC']] <- q$scDblFinder.score
CD28[['scDblFinder.class.ATAC']] <- q$scDblFinder.class

q <- scDblFinder(I$ATAC$counts, BPPARAM=bp)
I[['scDblFinder.score.ATAC']] <- q$scDblFinder.score
I[['scDblFinder.class.ATAC']] <- q$scDblFinder.class

CD28_Dblt_table = table(CD28[['scDblFinder.class.ATAC']]$scDblFinder.class.ATAC,CD28[['scDblFinder.class.RNA']]$scDblFinder.class.RNA)
print(CD28_Dblt_table)
fisher.test(CD28_Dblt_table)
mcnemar.test(CD28_Dblt_table)

I_Dblt_table = table(I[['scDblFinder.class.ATAC']]$scDblFinder.class.ATAC,I[['scDblFinder.class.RNA']]$scDblFinder.class.RNA)
print(I_Dblt_table)
fisher.test(I_Dblt_table)
mcnemar.test(I_Dblt_table)

# Create a list for Venn diagram
venn_list <- list(
  ATAC = rownames(CD28[['scDblFinder.class.ATAC']])[which(CD28[['scDblFinder.class.ATAC']]=="doublet")],
  RNA = rownames(CD28[['scDblFinder.class.RNA']])[which(CD28[['scDblFinder.class.RNA']]=="doublet")]
)
# Generate the Venn diagram
venn_plot <- ggvenn(venn_list, fill_color = c("blue", "red"))
ggsave(filename = "JOINT0_Venn_Doublet_CD28.pdf", plot = venn_plot, width = 6, height = 6)

# Create a list for Venn diagram
venn_list <- list(
  ATAC = rownames(I[['scDblFinder.class.ATAC']])[which(I[['scDblFinder.class.ATAC']]=="doublet")],
  RNA = rownames(I[['scDblFinder.class.RNA']])[which(I[['scDblFinder.class.RNA']]=="doublet")]
)
# Generate the Venn diagram
venn_plot <- ggvenn(venn_list, fill_color = c("blue", "red"))
ggsave(filename = "JOINT0_Venn_Doublet_Intact.pdf", plot = venn_plot, width = 6, height = 6)

CD28 <- subset(
  x = CD28,
  subset = scDblFinder.class.ATAC == "singlet" & scDblFinder.class.RNA == "singlet")

I <- subset(
  x = I,
  subset = scDblFinder.class.ATAC == "singlet" & scDblFinder.class.RNA == "singlet")

vln_plot <- VlnPlot(CD28, features = c("nFeature_RNA", "nCount_RNA","nCount_ATAC", "TSS.enrichment", "percent.mt"), ncol = 5, pt.size = 0.5)
ggsave(filename = "JOINT0_vln_plotCD28_Post.pdf", plot = vln_plot, width = 12, height = 6)

vln_plot <- VlnPlot(I, features = c("nFeature_RNA", "nCount_RNA","nCount_ATAC", "TSS.enrichment", "percent.mt"), ncol = 5, pt.size = 0.5)
ggsave(filename = "JOINT0_vln_plotI_Post.pdf", plot = vln_plot, width = 12, height = 6)


#UMAP

DefaultAssay(CD28) <- "RNA"
DefaultAssay(I) <- "RNA"

CD28 <- SCTransform(CD28)
CD28 <- RunPCA(CD28)
CD28 <- RunUMAP(CD28,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
pca_plot = DimPlot(CD28, reduction = "pca") + NoLegend()
ggsave(filename = "RNA1_pca_plotCD28.pdf", plot = pca_plot, width = 6, height = 6)
elbow_plot = ElbowPlot(CD28)
ggsave(filename = "RNA1_elbow_plotCD28.pdf", plot = elbow_plot, width = 6, height = 6)

I <- SCTransform(I)
I <- RunPCA(I)
I <- RunUMAP(I, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
pca_plot = DimPlot(I, reduction = "pca") + NoLegend()
ggsave(filename = "RNA1_pca_plotI.pdf", plot = pca_plot, width = 6, height = 6)
elbow_plot = ElbowPlot(I)
ggsave(filename = "RNA1_elbow_plotI.pdf", plot = elbow_plot, width = 6, height = 6)

DefaultAssay(CD28) <- "ATAC"
DefaultAssay(I) <- "ATAC"

CD28 <- FindTopFeatures(CD28, min.cutoff = 'q0')
CD28 <- RunTFIDF(CD28)
CD28 <- RunSVD(CD28)
CD28 = RunUMAP(CD28, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# build a joint neighbor graph using both assays
CD28 = FindMultiModalNeighbors(CD28, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
CD28 = RunUMAP(CD28, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CD28 = FindClusters(CD28, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

I <- FindTopFeatures(I, min.cutoff = 'q0')
I <- RunTFIDF(I)
I <- RunSVD(I)
I = RunUMAP(I, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# build a joint neighbor graph using both assays
I = FindMultiModalNeighbors(I, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
I = RunUMAP(I, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
I = FindClusters(I, graph.name = "wsnn", algorithm = 3, verbose = FALSE)



load(file = paste0("CD28_PP.RData"))
load(file = paste0("I_PP.RData"))

# Individual plots with titles and labels
p1 <- DimPlot(CD28, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle("RNA") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(CD28, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle("ATAC") + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(CD28, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle("WNN") + theme(plot.title = element_text(hjust = 0.5))

# Combine plots and remove legends
pfinal <- (p1 | p2 | p3) & NoLegend()

# Save combined plot to PDF
pdf("JOINT1_CD28_umaps.pdf", width = 15, height = 6)
print(pfinal)
dev.off()

# Individual plots with titles and labels
p1 <- DimPlot(I, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle("RNA") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(I, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle("ATAC") + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(I, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle("WNN") + theme(plot.title = element_text(hjust = 0.5))

# Combine plots and remove legends
pfinal <- (p1 | p2 | p3) & NoLegend()

# Save combined plot to PDF
pdf("JOINT1_I_umaps.pdf", width = 15, height = 6)
print(pfinal)
dev.off()

save(CD28, file = "CD28_PP.RData")
save(I, file = "I_PP.RData")


seu <- CD28 
save(seu, file = "CD28_PP.RData")

seu <- I 
save(seu, file = "I_PP.RData")





# Scoring / Cell typing

###################
## Author: Daniele Tavernari - daniele.tavernari@unil.ch
## Created: 18/04/2017
## Description: 
###################

####### START ######
rm(list = ls()) # clear all
#graphics.off() # close all
####################

#################### INPUT 

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

source("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/utils/dan.functions.R")

OutDir = "/mnt/etemp/ahrmad/wouter/signac/outDan"
dir.create(OutDir,showWarnings = F)

###############################################

################# Methods and Resources
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# https://stuartlab.org/signac/articles/pbmc_multiomic.html

# I am using a mixture of both resources.
# I will first process each sample separately, then maybe together (integrating?)
# For now I won't re-do peak calling (as suggested by Signac, but not by Seurat protocol). Maybe worth re-do peak calling on each cell subtype?

###############################################

#################### FUNCTIONS

###############################################

##### Each sample separately, for now
this_OutDir = paste0(OutDir,"/cell_typing/")
dir.create(this_OutDir)

# Sample = "WK-1585_INTACT_AP_BL6_Contrl" # 'AP' is Anterior, Proximal I guess. 'INTACT' and 'Control' means the initial states (before castration). 'BL6' is the type of mouse.
# # Also WK-1501_BL6_INTACT_AP_Test3 could be a good candidate to assess cell types

all_samples = c("I","CD28")

### Cell typing
# From Wouter presentation
# Cell types expected to be present: immune, mesenchymal, luminal, smooth muscle, basal (epi), neuroendocrine (epi), nerve
# Mouse: ratio epi:stromal = 50:50
# From Karthaus 2020 Science
# Luminal 1/2/3 (all markers in presentation/paper), Glial, Myeloid, Mesenchymal 1/2, Vascular endothelium, Lymphatic endothelium, T/B/NK/dendritic (few), myofibroblasts/smooth muscle
# Singulator protocol (data available right now) has different ratios of cell types. Much more L1 cells.
# From Methods section:
# splitting into epi (Epcam+), immune (Ptprc+), stromal (--). 
# B-cells, CD19, Ms4a1: T-cells, CD3, CD4, CD8; macrophages, CD14, Aif1; endothelium CD31, Vwf; lymphatic endothelium, CD31, Prox1; glia, Sox10; 
# myofibroblast, Acta2, Myh11, Rspo3; smooth muscle, Acta2, Notch3; mesenchymal 1, Col5a2, Lama2, Zeb1 Wnt2, Wnt6, Wnt10a, Rorb; 
# mesenschymal 2, Col5a2, Lama2, Zeb1, Sult1e1, Fgf10, Rspo1; basal, Epcam, Krt5, Krt14 Trp63; seminal vesicle basal, Epcam, Pax2, Krt5, Krt14 Trp63, Calml3; 
# luminal 1, Epcam, CD24a, Krt8, Krt18, Nkx3.1, Pbsn high; luminal 2, Epcam, CD24a, Krt8, Krt18, Psca, Krt4, Tacst2, Ly6a; luminal 3 (ionocyte), Epcam, CD24a, Krt8, Krt18, Foxi1, Atp6v1g3, Atp6b1b; 
# seminal vesicle luminal, Epcam, Pax2, Krt8, Krt18, Pate4; seminal vesicle ionocyte, Epcam, Pax2, Foxi1.

main_ct_marker_list = list(Epithelial = c("Epcam"), Immune = "Ptprc") # present
# sub_ct_marker_list = list(Bcells = c("Cd19","Ms4a1"), Tcells = c("Cd3g","Cd3d","Cd3e","Cd247"), macrophages = c("Cd14", "Aif1"), endothelium = c("Pecam1","Vwf"),
#       lymphatic_endothelium = c("Pecam1","Prox1" ), glia = c("Sox10"), myofibroblast = c("Acta2","Myh11","Rspo3" ), smooth_muscle = c("Acta2", "Notch3"),
#       mesenchymal_1 = c("Col5a2", "Lama2", "Zeb1", "Wnt2", "Wnt6", "Wnt10a", "Rorb"), mesenschymal_2 = c("Col5a2", "Lama2", "Zeb1", "Sult1e1", "Fgf10", "Rspo1"),
#       basal = c("Epcam", "Krt5", "Krt14", "Trp63"), seminal_vesicle_basal = c("Epcam", "Pax2", "Krt5", "Krt14","Trp63", "Calml3"), 
#       luminal_1 = c("Epcam", "Cd24a", "Krt8", "Krt18", "Nkx3-1", "Pbsn", PROM1, CD26/DPP4), luminal_2 = c("Epcam", "Cd24a", "Krt8", "Krt18", "Psca", "Krt4", "Tacstd2", "Ly6a"),
#       luminal_3 = c("Epcam", "Cd24a", "Krt8", "Krt18", "Foxi1", "Atp6v1g3", "Atp6v1b1"), seminal_vesicle_luminal = c( "Epcam", "Pax2", "Krt8", "Krt18", "Pate4"),
#        seminal_vesicle_ionocyte = c("Epcam", "Pax2", "Foxi1") ) # already alias-corrected
### big absences from 'measured_genes_all': none
### big absences from 'measured_genes' (sample-specific): Cd19, Cd3*, Cd31/Pecam1, Sox10, some Wnt, Cd24a, ("Foxi1", "Atp6v1g3", "Atp6v1b1") discriminating the luminal_3, 

sub_ct_marker_list_table_s8 = dan.read("/mnt/ndata/daniele/wouter/Data/table_s8_summary.txt")
sub_ct_marker_list_table_s8 = as.list(sub_ct_marker_list_table_s8)
sub_ct_marker_list_table_s8 = lapply(sub_ct_marker_list_table_s8, function(x) x = x[nchar(x)>0])
print(names(sub_ct_marker_list_table_s8))

sub_ct_marker_list_table_s9 = dan.read("/mnt/ndata/daniele/wouter/Data/table_s9_summary.txt")
sub_ct_marker_list_table_s9 = as.list(sub_ct_marker_list_table_s9)
sub_ct_marker_list_table_s9 = lapply(sub_ct_marker_list_table_s9, function(x) x = x[nchar(x)>0])
print(names(sub_ct_marker_list_table_s9))

# j89 = dan.df(names(sub_ct_marker_list_table_s8),names(sub_ct_marker_list_table_s9))
# for (ii in names(sub_ct_marker_list_table_s8)){
#   for (jj in names(sub_ct_marker_list_table_s9)){
#     j89[ii,jj] = length(intersect(sub_ct_marker_list_table_s8[[ii]],sub_ct_marker_list_table_s9[[jj]]))/length(union(sub_ct_marker_list_table_s8[[ii]],sub_ct_marker_list_table_s9[[jj]]))
#   }
# }


#Sample = "WK-1585_INTACT_AP_BL6_Contrl" # 'AP' is Anterior, Proximal I guess. 'INTACT' and 'Control' means the initial states (before castration). 'BL6' is the type of mouse.
# Also WK-1501_BL6_INTACT_AP_Test3 could be a good candidate to assess cell types

for (Sample in all_samples){

  dcat(Sample)

  load(file = paste0(Sample,"_PP.RData"))
  measured_genes_all = rownames(seu@assays$RNA)
  measured_genes = rownames(seu@assays$SCT) 

  # lapply(sub_ct_marker_list, function(x) intersect(x, measured_genes_all) )
  # lapply(sub_ct_marker_list, function(x) intersect(x, measured_genes) )
  # "Atp6v1b1" %in% measured_genes_all
  # "Cd24a" %in% measured_genes_all

  ## Main cell type markers
  DefaultAssay(seu) = "SCT"
  seu = FindNeighbors(seu)
  seu = FindClusters(seu,res=1.5)
  these_features = intersect(as.character(unlist(main_ct_marker_list)),rownames(seu))
  p1 = FeaturePlot(seu, features = these_features, reduction = "umap.rna",  label = TRUE, label.size = 2.5, ncol = length(these_features) ) + coord_fixed()
  pdf( paste0(this_OutDir,Sample,"_MainCellTypes_umaps.pdf"),length(these_features)*4,4 )
  patchwork = wrap_plots(p1, nrow = 1)
  print(patchwork)
  dev.off()

  ### Signature scoring, table s8
  these_sub_ct_marker_list = sub_ct_marker_list_table_s8[names(which(sapply(sub_ct_marker_list_table_s8, function(x) length(intersect(x, measured_genes)) )>0))]
  thisseu = AddModuleScore(seu, features = these_sub_ct_marker_list, name = "subct")
  colnames(thisseu@meta.data)[substr(colnames(thisseu@meta.data),1,5)=="subct"] = names(these_sub_ct_marker_list)
  these_features = names(these_sub_ct_marker_list)
  p1 = FeaturePlot(thisseu, features = these_features, reduction = "umap.rna",  label = TRUE, label.size = 2.5, ncol = 7) + coord_fixed()
  pdf( paste0(this_OutDir,Sample,"_SubCellTypeScores_TableS8_umaps.pdf"),24,9 )
  print(p1)
  dev.off()

  colorz = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))
  pdf( paste0(this_OutDir,Sample,"_SubCellTypeScores_TableS8_umaps_SameScale.pdf"),24,9 )
  p1 = FeaturePlot(thisseu, features = these_features, reduction = "umap.rna",  label = TRUE, label.size = 2.5, ncol = 7, keep.scale = "all", cols = colorz) + coord_fixed()
  print(p1)
  dev.off()

  these_features = c("Epi_Basal_1","Epi_Basal_SV","Epi_Luminal_1","Epi_Luminal_2Psca","Epi_Luminal_3Foxi1","Epi_Luminal_SV"  )
  colorz = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))
  pdf( paste0(this_OutDir,Sample,"_SubCellTypeScores_TableS8_umaps_SameScale_Epi.pdf"),15,9 )
  p1 = FeaturePlot(thisseu, features = these_features, reduction = "umap.rna",  label = TRUE, label.size = 2.5, ncol = 3, keep.scale = "all", cols = colorz) + coord_fixed()
  print(p1)
  dev.off()

  # ### Signature scoring, table s9
  # these_sub_ct_marker_list = sub_ct_marker_list_table_s9[names(which(sapply(sub_ct_marker_list_table_s9, function(x) length(intersect(x, measured_genes)) )>0))]
  # thisseu = AddModuleScore(seu, features = these_sub_ct_marker_list, name = "subct")
  # colnames(thisseu@meta.data)[substr(colnames(thisseu@meta.data),1,5)=="subct"] = names(these_sub_ct_marker_list)
  # these_features = names(these_sub_ct_marker_list)
  # p1 = FeaturePlot(thisseu, features = these_features, reduction = "umap.rna",  label = TRUE, label.size = 2.5, ncol = 6) + coord_fixed()
  # pdf( paste0(this_OutDir,Sample,"_SubCellTypeScores_TableS9_umaps.pdf"),21,9 )
  # print(p1)
  # dev.off()

  # colorz = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))
  # pdf( paste0(this_OutDir,Sample,"_SubCellTypeScores_TableS9_umaps_SameScale.pdf"),21,9 )
  # p1 = FeaturePlot(thisseu, features = these_features, reduction = "umap.rna",  label = TRUE, label.size = 2.5, ncol = 7, keep.scale = "all", cols = colorz) + coord_fixed()
  # print(p1)
  # dev.off()

  md = seu@meta.data
  save(md, file = paste0(this_OutDir,Sample,"_clusters.RData"))

}

