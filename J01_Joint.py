import matplotlib.pyplot as plt
import snapatac2 as snap
import scanpy as sc
import pandas as pd
import numpy as np
import os

for met in ['SING','ENZ']:
	refDir = '/mnt/etemp/ahrmad/wouter/refs'
	### TO CHECK WHEN CHANGING SAMPLES ###
	DirRNA = f'/mnt/etemp/ahrmad/wouter/{met}_RNA_CB'
	DirATAC = f'/mnt/etemp/ahrmad/wouter/{met}_ATAC'
	Exp_ATAC = f'Wouter21_{met}'
	Exp_RNA = f'Wouter21_{met}_CB'
	Experiment=f'Wouter21_{met}'
	### TO CHECK WHEN CHANGING SAMPLES ###
	gene_set_table = f'{refDir}/table_s8_summary.txt'


	if not (os.path.exists(f'{Experiment}_GEX_CB_GeneScores.h5ad')) & (os.path.exists(f'{Experiment}_ATAC_GeneScores.h5ad')):
		print('Making gene scored H5ad files...')
		#rna = snap.read(DirRNA +'/'+ Exp_RNA + '_post.h5ad').to_memory() #notimplemented
		#atac = snap.read(DirATAC +'/'+ Exp_ATAC + '_Post.h5ad').to_memory()

		rna = sc.read_h5ad(DirRNA +'/'+ Exp_RNA + '_post.h5ad')
		#rna = sc.read_h5ad(DirRNA +'/'+ Exp_RNA + '_Norm_FeatSel_DimRed.h5ad')
		atac = sc.read_h5ad(DirATAC +'/'+ Exp_ATAC + '_GeneMat.h5ad')
		
		#Celltyping genes based on Supp gene lists
		gset_df = pd.read_csv(gene_set_table, sep='\t')
		gset_dict = dict()
		for gset_name in gset_df.columns:
			gset_dict[gset_name] = gset_df[gset_name].dropna().str.upper().tolist()
		gset_dict = {key: [item for item in value if item in common_genes] for key, value in gset_dict.items()}

		# Make sure the same cells are in both adatas
		print(f'keeping cells in both modalities')
		atac.obs_names = [el.split('_')[-1]+'_'+'_'.join(el.split('_')[0:-1]) for el in atac.obs_names]
		rna.obs = rna.obs.reindex(index=sorted(rna.obs.index))
		atac.obs = atac.obs.reindex(index=sorted(atac.obs.index))
		lenATAC,lenGEX = atac.n_obs,rna.n_obs
		atac = atac[atac.obs.index.isin(rna.obs.index)].copy()
		rna = rna[rna.obs.index.isin(atac.obs.index)].copy()
		assert (rna.obs.index == atac.obs.index).all()
		print(f'rna lost {lenGEX-rna.n_obs} cells\natac lost {lenATAC-atac.n_obs} cells')
		
		# only keep genes that are present in both ATAC and GEX
		print(f'finding genes present in both modalities')
		lenGEX,lenATAC = rna.n_vars,atac.n_vars
		rna.var.index = rna.var.index.str.upper()
		atac.var.index = atac.var.index.str.upper()
		common_genes = rna.var.index.intersection(atac.var.index).tolist()
		print(f'rna lost {lenGEX-len(common_genes)} genes\natac lost {lenATAC-len(common_genes)} genes')

		# Cell gene scoring
		for gset_name in gset_dict.keys():
			gs = gset_dict[gset_name]
			print(f'{gset_name}: {len(gs)}genes Karthaus2020 TableS8')
			
			ctrl = len(gs) if (len(gs) > 50) else 50 #genes
			sc.tl.score_genes(rna, gs, ctrl_size=ctrl, n_bins=25, score_name=gset_name,use_raw=False)
			sc.tl.score_genes(atac, gs, ctrl_size=ctrl, n_bins=25, score_name=gset_name,use_raw=False)

		rna.write(f'{Experiment}_GEX_CB_GeneScores.h5ad')
		atac.write(f'{Experiment}_ATAC_GeneScores.h5ad')
	
	else:
		print('Reading H5ad files...')
		rna = sc.read_h5ad(f'{Experiment}_GEX_CB_GeneScores.h5ad')
		atac = sc.read_h5ad(f'{Experiment}_ATAC_GeneScores.h5ad')
		# only keep genes that are present in both ATAC and GEX
		print(f'finding genes present in both modalities')
		lenGEX,lenATAC = rna.n_vars,atac.n_vars
		rna.var.index = rna.var.index.str.upper()
		atac.var.index = atac.var.index.str.upper()
		common_genes = rna.var.index.intersection(atac.var.index).tolist()
		print(f'rna lost {lenGEX-len(common_genes)} genes\natac lost {lenATAC-len(common_genes)} genes')

		#Celltyping genes based on Supp gene lists
		gset_df = pd.read_csv(gene_set_table, sep='\t')
		gset_dict = dict()
		for gset_name in gset_df.columns:
			gset_dict[gset_name] = gset_df[gset_name].dropna().str.upper().tolist()
		gset_dict = {key: [item for item in value if item in common_genes] for key, value in gset_dict.items()}


	minATAC,maxATAC,minRNA,maxRNA = [],[],[],[]
	for gset_name in gset_dict.keys():
		minATAC.append(np.min(atac.obs[gset_name]))
		maxATAC.append(np.max(atac.obs[gset_name]))
		minRNA.append(np.min(rna.obs[gset_name]))
		maxRNA.append(np.max(rna.obs[gset_name]))
	# Round by 0.5 below or above
	#atacMinMax = np.floor(min(minATAC) * 2) / 2, np.ceil(max(maxATAC) * 2) / 2 
	#rnaMinMax = np.floor(min(minRNA) * 2) / 2, np.ceil(max(maxRNA) * 2) / 2
	atacMinMax = min(minATAC)-0.1,max(maxATAC)+0.1
	rnaMinMax = min(minRNA)-0.1,max(maxRNA)+0.1

	def scatter_gene_sigs(rna,atac,gsetName,Experiment=Experiment,atacylim=atacMinMax,rnaxlim=rnaMinMax):
		assert (rna.obs.index == atac.obs.index).all()
		# Create a scatter plot with colors based on the average sig score
		sc = plt.scatter(rna.obs[gsetName], atac.obs[gsetName],
		c=(rna.obs[gsetName] + atac.obs[gsetName])/2, cmap='viridis', s=6)

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


	for gset_name in gset_dict.keys():
		print(f'{Experiment}: {gset_name} scatter plot')
		scatter_gene_sigs(rna,atac,gset_name)

#pdfjam --nup 4x2 --landscape --a4paper *ENZ*average*Imm*.png --outfile ENZ_Imm_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper *ENZ*average*Str*.png --outfile ENZ_Str_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper *ENZ*average*Epi*.png --outfile ENZ_Epi_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper *SING*average*Imm*.png --outfile SING_Imm_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper *SING*average*Str*.png --outfile SING_Str_CB.pdf
#pdfjam --nup 4x2 --landscape --a4paper *SING*average*Epi*.png --outfile SING_Epi_CB.pdf





