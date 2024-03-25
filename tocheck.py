import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import snapatac2 as snap
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import os

#for met in ['SING','ENZ']:

met='SING'
refDir = '/mnt/etemp/ahrmad/wouter/refs'
### TO CHECK WHEN CHANGING SAMPLES ###
DirRNA = f'/mnt/etemp/ahrmad/wouter/{met}_RNA_CB'
DirATAC = f'/mnt/etemp/ahrmad/wouter/{met}_ATAC'
Exp_ATAC = f'Wouter21_{met}'
Exp_RNA = f'Wouter21_{met}_CB'
Experiment=f'Wouter21_{met}'
### TO CHECK WHEN CHANGING SAMPLES ###
annot_ext = '_Annotated.h5ad'
gene_set_table = f'{refDir}/table_s8_summary.txt'

print('Reading H5ad files...')
rna = sc.read_h5ad(f'{Experiment}_RNA{annot_ext}')
atac = sc.read_h5ad(f'{Experiment}_ATAC{annot_ext}')

resLeiden = 0.1
sc.tl.leiden(rna, key_added=f"leiden_res{resLeiden}", resolution=resLeiden)
#sc.tl.leiden(atac, key_added=f"leiden_res{resLeiden}", resolution=resLeiden)

print(f'keeping cells in both modalities')
atac.obs_names = [el.split('_')[-1]+'_'+'_'.join(el.split('_')[0:-1]) for el in atac.obs_names]
rna.obs = rna.obs.reindex(index=sorted(rna.obs.index))
atac.obs = atac.obs.reindex(index=sorted(atac.obs.index))
lenATAC,lenGEX = atac.n_obs,rna.n_obs
atac = atac[atac.obs.index.isin(rna.obs.index)].copy()
rna = rna[rna.obs.index.isin(atac.obs.index)].copy()
assert (rna.obs.index == atac.obs.index).all()
print(f'rna lost {lenGEX-rna.n_obs} cells\natac lost {lenATAC-atac.n_obs} cells')

#Celltyping genes based on Supp gene lists
gset_df = pd.read_csv(gene_set_table, sep='\t')
gset_dict = dict()
for gset_name in gset_df.columns:
	gset_dict[gset_name] = gset_df[gset_name].dropna().str.upper().tolist()



def scatter_gene_sigs(rna,atac,gsetName,Experiment=Experiment):
	assert (rna.obs.index == atac.obs.index).all()

	minATAC,maxATAC,minRNA,maxRNA = [],[],[],[]
	for gset_name in gset_dict.keys():
		minATAC.append(np.min(atac.obs[gset_name]))
		maxATAC.append(np.max(atac.obs[gset_name]))
		minRNA.append(np.min(rna.obs[gset_name]))
		maxRNA.append(np.max(rna.obs[gset_name]))
	atacylim = min(minATAC)+0.04,max(maxATAC)-0.2
	rnaxlim = min(minRNA)+0.2,max(maxRNA)-0.2

	# Create a scatter plot with colors based on the average sig score
	sc = plt.scatter(rna.obs[gsetName], atac.obs[gsetName],
	c=(rna.obs[gsetName] + atac.obs[gsetName])/2, cmap='viridis', s=0.2)

	# Label axes and set a title
	plt.xlabel(f'GEX {gsetName} Score')
	plt.ylabel(f'ATAC {gsetName} Score')
	plt.title(f'{gsetName} Score in ATAC and GEX')
	plt.suptitle(Experiment)
	
	# Display the Pearson correlation coefficient in a box in the top right corner
	plt.text(0.85, 0.9, f'r = {rna.obs[gsetName].corr(atac.obs[gsetName]):.2f}', transform=plt.gca().transAxes,
		bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='round,pad=0.3'), fontsize=10)

	# Set x and y-axis limits
	plt.xlim(rnaxlim[0], rnaxlim[1])  # Adjust the limits based on data
	plt.ylim(atacylim[0], atacylim[1])  # Adjust the limits based on data

	# Display the plot
	plt.tight_layout()

	# Save the scatter plot
	plt.savefig(f'{Experiment}_average_{gsetName}_Score.png')
	plt.close()

def boxplot_gene_sigs(mod,gsetName,modality,Experiment=Experiment,clusters='leiden_res0.25'):
	
	# Create a boxplot using seaborn
	sns.set(style="whitegrid")

	plt.figure(figsize=(10, 6))
	sns.boxplot(x=mod.obs[clusters], y=mod.obs[gsetName], hue=mod.obs[clusters], dodge=False)

	# Add titles and labels
	plt.title(f'Boxplot of {gsetName} Scores by {clusters} - {Experiment}')
	plt.xlabel(f'Clusters {clusters}')
	plt.ylabel(f'{gsetName} Score')

	plt.tight_layout()
	plt.savefig(f'{Experiment}_{modality}boxplots_{clusters}_{gsetName}_Score.png')
	plt.close()


for gset_name in ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Luminal_3Foxi1']:
	print(f'{Experiment}: {gset_name} scatter plot')
	scatter_gene_sigs(rna,atac,gset_name)

	print(f'{Experiment}: {gset_name} boxplots')
	boxplot_gene_sigs(atac,gset_name,'ATAC',clusters='leiden_res0.25')
	boxplot_gene_sigs(rna,gset_name,'RNA',clusters='leiden_res0.1')

rna_cl7 = rna[rna.obs['leiden_res0.25'] == '7']
atac_cl7 = atac[atac.obs.index.isin(rna_cl7.obs.index)].copy()
assert (rna_cl7.obs.index == atac_cl7.obs.index).all()

for gset_name in ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca','Epi_Luminal_3Foxi1']:
	scatter_gene_sigs(rna_cl7,atac_cl7,gset_name,Experiment=f'Wouter21_{met}_cl7')

rna_cl0 = rna[rna.obs['leiden_res0.25'] == '0']
atac_cl0 = atac[atac.obs.index.isin(rna_cl0.obs.index)].copy()
assert (rna_cl0.obs.index == atac_cl0.obs.index).all()

for gset_name in gset_df.columns:
	print(f'{Experiment}: {gset_name} scatter plot')
	scatter_gene_sigs(rna_cl0,atac_cl0,gset_name,Experiment=f'Wouter21_{met}_cl0')

pdfjam --nup 4x1 --landscape --a4paper *ENZ*average*Epi*.png --outfile ENZ_Epi_CB.pdf
pdfjam --nup 4x1 --landscape --a4paper *SING*average*Epi*.png --outfile SING_Epi_CB.pdf

pdfjam --nup 2x2 --landscape --a4paper *ENZ*ATACboxplots*Epi*.png --outfile ENZ_Epi_ATAC_box.pdf
pdfjam --nup 2x2 --landscape --a4paper *ENZ*RNAboxplots*Epi*.png --outfile ENZ_Epi_RNA_box.pdf

pdfjam --nup 2x2 --landscape --a4paper *SING*ATACboxplots*Epi*.png --outfile SING_Epi_ATAC_box.pdf
pdfjam --nup 2x2 --landscape --a4paper *SING*RNAboxplots*Epi*.png --outfile SING_Epi_RNA_box.pdf

#SPECIAL
pdfjam --nup 4x1 --landscape --a4paper *ENZ_cl7*average*Epi*.png --outfile ENZ_cl7_Epi_CB.pdf
pdfjam --nup 7x3 --landscape --a4paper *ENZ_cl0*average*.png --outfile ENZ_cl0_Epi_CB.pdf

pdfjam --nup 2x2 --landscape --a4paper *SING*RNAboxplots_leiden_res0.1*Epi*.png --outfile SING_Epi_RNA_box_res0_1.pdf



import snapatac2 as snap

atac_re = atac.copy()
snap.pp.select_features(atac_re, n_features=3000, inplace=True) #

#Perform dimension reduction using the spectrum of the normalized graph Laplacian defined by pairwise similarity between cells
snap.tl.spectral(atac_re, weighted_by_sd=True, features='selected', distance_metric='cosine', inplace=True)
snap.tl.umap(atac_re, use_rep='X_spectral', key_added='umap', random_state=None)

#neighborhood graph of observations stored in data using the method specified by method. The distance metric used is Euclidean.
snap.pp.knn(atac_re, n_neighbors=50, use_rep='X_spectral', method='kdtree')

#Cluster cells using the Leiden algorithm [Traag18]
for resLeiden in [.25,.5,1,1.5,2]:
	print(f'Leiden clustering at {resLeiden} resolution')
	snap.tl.leiden(atac_re, resolution=resLeiden, key_added=f"leiden_res{resLeiden}")
	snap.pl.umap(atac_re, color=f"leiden_res{resLeiden}", height=500,interactive=False, show=False, out_file=f"{Experiment}_leiden_res{resLeiden}.pdf")
	sc.pl.umap(atac_re,color=f"leiden_res{resLeiden}",legend_loc="on data",save=f"{Experiment}_leiden_res{resLeiden}_sc.pdf",title=Experiment,show=False)






