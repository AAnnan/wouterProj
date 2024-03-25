import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import snapatac2 as snap
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
	annot_ext = '_Annotated.h5ad'
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

		rna.write(Experiment+annot_ext)
		atac.write(Experiment+annot_ext)
	
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


### Scatter Gene Signatures
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

### Cluster analysis

rna_leiden_Res='leiden_res0.25'
atac_leiden_Res='leiden_res0.25'

for met in ['SING','ENZ']:
	Experiment=f'Wouter21_{met}'
	print(f'Reading {Experiment} H5ad files...')
	rna = sc.read_h5ad(f'{Experiment}_RNA{annot_ext}')
	atac = sc.read_h5ad(f'{Experiment}_ATAC{annot_ext}')

	for sample in rna.obs['batch'].cat.categories:
		dataGEX = rna[rna.obs['batch']==sample].copy()
		dataATAC = atac[atac.obs['batch']==sample].copy()

		GEXdf = pd.DataFrame({'sample': sample, 'obs_names': dataGEX.obs_names, 'leiden': dataGEX.obs[rna_leiden_Res]})
		ATACdf = pd.DataFrame({'sample': sample, 'obs_names': dataATAC.obs_names, 'leiden': dataATAC.obs[atac_leiden_Res]})

		merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')
		#commonBClist = list(set(dataATAC.obs.index).intersection(set(dataGEX.obs.index)))
		#print(np.all(np.array([True for el in commonBClist if el in list(merged_df['obs_names'])])),'Barcode Intersection')

		print(f'Plotting {sample}')
		#########
		#Matrix comparison
		#########

		groups1 = np.array(merged_df['leiden_x'],dtype=np.int16)
		cl1 = np.unique(groups1)
		groups2 = np.array(merged_df['leiden_y'],dtype=np.int16)
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
		fig, ax = plt.subplots()
		im = ax.imshow(heatmap_data, cmap='coolwarm', vmin=0, vmax=1)
	    
	    # set up the plot
		ax.set_xticks(np.arange(len(cl2)))
		ax.set_yticks(np.arange(len(cl1)))
		#Put ticks in right order 
		ax.set_xticklabels(cl2[sorted_indices_x])
		ax.set_yticklabels(cl1[sorted_indices_y])

	    # Add the number of observations in each group on x and y axes
	    for i, (value, label) in enumerate(zip(cl2[sorted_indices_x], ax.get_xticklabels())):
	        ax.text(i, -0.56, f'{len(np.where(groups2 == value)[0])}', ha='left', va='center',rotation=45,rotation_mode='anchor')
	    
	    for i, (value, label) in enumerate(zip(cl1[sorted_indices_y], ax.get_yticklabels())):
	        ax.text(-0.45, i, f'{len(np.where(groups1 == value)[0])}', ha='left', va='center')
	        label.set_y(i - 0.3)  # Adjust the position of the current label

		ax.set_xlabel('GEX Clusters')
		ax.set_ylabel('ATAC Clusters')
		ax.set_title(f'{sample}',pad=30)
		divider = make_axes_locatable(plt.gca())
		cax = divider.append_axes("right", "5%", pad="3%")
		plt.colorbar(im, cax=cax)
		plt.tight_layout()

		# save the plot
		aa = atac_leiden_Res.split('leiden_res')[1]
		rr = rna_leiden_Res.split('leiden_res')[1]
		plt.savefig(f'{Experiment}_{sample}_a{aa}_r{rr}_Jaccard.pdf')
		plt.close()

pdfjam --nup 4x3 --landscape --a4paper *ENZ*WK*a0.25*Jaccard.pdf --outfile ENZ_Jacc_a0.25.pdf
pdfjam --nup 4x3 --landscape --a4paper *SING*WK*a0.25*Jaccard.pdf --outfile SING_Jacc_a0.25.pdf

### Cluster analysis pan samples
rna_leiden_Res='leiden_res0.1'
atac_leiden_Res='leiden_res0.25'
for met in ['SING','ENZ']:
	Experiment=f'Wouter21_{met}'
	print(f'Reading {Experiment} H5ad files...')
	rna = sc.read_h5ad(f'{Experiment}_GEX_CB_GeneScores.h5ad')
	atac = sc.read_h5ad(f'{Experiment}_ATAC_GeneScores.h5ad')

	GEXdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': rna.obs_names, 'leiden': rna.obs[rna_leiden_Res]})
	ATACdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': atac.obs_names, 'leiden': atac.obs[atac_leiden_Res]})
	merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')
	#commonBClist = list(set(atac.obs.index).intersection(set(rna.obs.index)))
	#print(np.all(np.array([True for el in commonBClist if el in list(merged_df['obs_names'])])),'Barcode Intersection')

	print(f'Plotting {Experiment}')
	#########
	#Matrix comparison
	#########

	groups1 = np.array(merged_df['leiden_x'],dtype=np.int16)
	cl1 = np.unique(groups1)
	groups2 = np.array(merged_df['leiden_y'],dtype=np.int16)
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
	fig, ax = plt.subplots()
	im = ax.imshow(heatmap_data, cmap='coolwarm', vmin=0, vmax=1)

	# set up the plot
	ax.set_xticks(np.arange(len(cl2)))
	ax.set_yticks(np.arange(len(cl1)))
	#Put ticks in right order 
	ax.set_xticklabels(cl2[sorted_indices_x])
	ax.set_yticklabels(cl1[sorted_indices_y])

	# Add the number of observations in each group on x and y axes
	for i, (value, label) in enumerate(zip(cl2[sorted_indices_x], ax.get_xticklabels())):
	    ax.text(i, -0.56, f'{len(np.where(groups2 == value)[0])}', ha='left', va='center',rotation=45,rotation_mode='anchor')

	for i, (value, label) in enumerate(zip(cl1[sorted_indices_y], ax.get_yticklabels())):
	    ax.text(-0.45, i, f'{len(np.where(groups1 == value)[0])}', ha='left', va='center')
	    label.set_y(i - 0.3)  # Adjust the position of the current label

	ax.set_xlabel('GEX Clusters')
	ax.set_ylabel('ATAC Clusters')
	ax.set_title(f'{Experiment}',pad=30)
	divider = make_axes_locatable(plt.gca())
	cax = divider.append_axes("right", "5%", pad="3%")
	plt.colorbar(im, cax=cax)
	plt.tight_layout()

	# save the plot
	aa = atac_leiden_Res.split('leiden_res')[1]
	rr = rna_leiden_Res.split('leiden_res')[1]
	plt.savefig(f'{Experiment}_pan_a{aa}_r{rr}_Jaccard.pdf')
	plt.close()


# L1/L2 Signature
def karthaus_plotting(adata,timepoint,Experimenti,res_leid,cL1,cL2):
    x = adata.obs['Epi_Luminal_1']
    y = adata.obs['Epi_Luminal_2Psca']
    colLuminal = []
    for val in adata.obs[res_leid]:
        if val in cL1:
            colLuminal.append('blue')
        elif val in cL2:
            colLuminal.append('red')
        else:
            colLuminal.append('gray') 

    # Create a scatter plot with different colors for x and y values
    plt.scatter(x, y,c=colLuminal, s=8,alpha=0.5)
    # Add labels and title
    plt.xlabel('L1 Signature Score')
    plt.ylabel('L2 Signature Score')
    plt.title(f'{Experimenti} {timepoint}')
    
    #legend_labels = {label: f'L{i}' for i, labels in enumerate([cL1, cL2], 1) for label in labels}
    #handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=label)
    #           for color, label in legend_labels.items()]
    #plt.legend(handles=handles, title='Legend', loc='upper right')

    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"L1L2Sig_{Experimenti}_{timepoint}2.pdf")
    plt.close()
    print(f"Plotted L1/L2 {Experimenti} {timepoint}")

Experiment = 'Wouter21_ENZ'

adata = sc.read_h5ad(f'{Experiment}_GEX_CB_GeneScores.h5ad')
i_adata = adata[adata.obs['batch'].str.contains('WK-1350_I')].copy()
r3_adata = adata[adata.obs['batch'].str.contains('WK-1350_R3')].copy()

res_leiden = 'leiden_res1'
cl_L1 = ['7']
cl_L2 = ['16']

karthaus_plotting(i_adata,'Intact',Experiment,res_leiden,cl_L1,cl_L2)
karthaus_plotting(r3_adata,'RegenDay3',Experiment,res_leiden,cl_L1,cl_L2)

Experiment = 'Wouter21_SING'
res_leiden = 'leiden_res1'
cl_L1 = ['22','6','16','2']
cl_L2 = ['19']

adata = sc.read_h5ad(f'{Experiment}_GEX_CB_GeneScores.h5ad')

i_adata = adata[adata.obs['batch'].str.upper().str.contains('INTACT')].copy()
r3_adata = adata[adata.obs['batch'].str.contains('Day3')].copy()
r2_adata = adata[adata.obs['batch'].str.contains('Day2')].copy()
r1_adata = adata[adata.obs['batch'].str.contains('Day1')].copy()
c28_adata = adata[adata.obs['batch'].str.contains('Day28')].copy()

karthaus_plotting(i_adata,'Intact',Experiment,res_leiden,cl_L1,cl_L2)
karthaus_plotting(r3_adata,'RegenDay3',Experiment,res_leiden,cl_L1,cl_L2)
karthaus_plotting(r2_adata,'RegenDay2',Experiment,res_leiden,cl_L1,cl_L2)
karthaus_plotting(r1_adata,'RegenDay1',Experiment,res_leiden,cl_L1,cl_L2)
karthaus_plotting(c28_adata,'CastrateDay28',Experiment,res_leiden,cl_L1,cl_L2)







