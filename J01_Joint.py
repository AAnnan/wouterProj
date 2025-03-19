import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import snapatac2 as snap
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde
import os

### TO CHECK WHEN CHANGING SAMPLES ###
met = 'SING'#,'SING']:
### TO CHECK WHEN CHANGING SAMPLES ###


refDir = '/mnt/etemp/ahrmad/wouter/refs'
DirRNA = f'/mnt/etemp/ahrmad/wouter/RNA_CB/annot_other_norm/'
DirATAC = f'/mnt/etemp/ahrmad/wouter/Wouter21_{met}/annot/'
Exp_ATAC = f'Wouter21_{met}'
Exp_RNA = f'Wouter21_{met}_CB'
Experiment=f'Wouter21_{met}'
annot_ext = '_annot_Lum.h5ad'
gene_set_table = f'{refDir}/table_s8_summary.txt'
gene_sigs = ['Epi_Luminal_1', 'Epi_Luminal_2Psca']
NoSORT = True

assert (os.path.exists(DirRNA+ Exp_RNA + annot_ext)) & (os.path.exists(DirATAC+ Exp_ATAC + annot_ext)),'ScoreGene H5ADs missing'

rna = sc.read_h5ad(DirRNA+ Exp_RNA + annot_ext)
atac = sc.read_h5ad(DirATAC+ Exp_ATAC + annot_ext)

if 'SING' in Experiment and NoSORT:
    rna = rna[~(rna.obs['batch']=='WK-1501_BL6_INTACT_AP_Test3_SORT')].copy()
    atac = atac[~(atac.obs['batch']=='WK-1501_BL6_INTACT_AP_Test3_SORT')].copy()

# Make sure the same cells are in both adatas
print(f'keeping cells in both modalities')
rna.obs = rna.obs.reindex(index=sorted(rna.obs.index))
atac.obs = atac.obs.reindex(index=sorted(atac.obs.index))
lenATAC,lenGEX = atac.n_obs,rna.n_obs
atac = atac[atac.obs.index.isin(rna.obs.index)]
rna = rna[rna.obs.index.isin(atac.obs.index)]
assert (rna.obs.index == atac.obs.index).all()
print(f'rna lost {lenGEX-rna.n_obs} cells\natac lost {lenATAC-atac.n_obs} cells')

### Scatter/Density Gene Signatures

def scatter_gene_sigs(rna,atac,gsetName,dens,Experiment=Experiment):
    assert (rna.obs.index == atac.obs.index).all()
    if dens:
	    # Calculate the kernel density estimate
	    kde = gaussian_kde([rna.obs[gsetName], atac.obs[gsetName]], bw_method=0.1)

	    # Create a grid of points for the density estimation
	    x, y = np.mgrid[rna.obs[gsetName].min():rna.obs[gsetName].max():100j, atac.obs[gsetName].min():atac.obs[gsetName].max():100j]
	    positions = np.vstack([x.ravel(), y.ravel()])
	    t = np.reshape(kde(positions).T, x.shape)
	    ##z = np.log2(t)
	    z = np.sqrt(t)
	    xmin,xmax = np.min(x),np.max(x)    
    	ymin,ymax = np.min(y),np.max(y)
    else:
    	xmin,xmax = np.min(rna.obs[gsetName]),np.max(rna.obs[gsetName])    
    	ymin,ymax = np.min(atac.obs[gsetName]),np.max(atac.obs[gsetName])

    # Create a grid for the figure layout
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(3, 3, width_ratios=[6, 1, 0.1], height_ratios=[0.1, 6, 1])

    # Main scatter plot
    main_ax = fig.add_subplot(gs[1, 0])
    if dens:
    	main_ax.contourf(x, y, z, cmap='plasma', levels=40, antialiased=True)  # Create a filled contour plot for density
    else:
    	main_ax.scatter(rna.obs[gsetName], atac.obs[gsetName], c='#377eb8', marker='o', s=10, edgecolors='#377eb8')  # Scatter plot

    main_ax.text(0.9, 0.95, f'r = {rna.obs[gsetName].corr(atac.obs[gsetName]):.2f}', transform=main_ax.transAxes,
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='round,pad=0.3'), fontsize=12)
    main_ax.set_xlabel(f'GEX {gsetName} Score')
    main_ax.set_ylabel(f'ATAC {gsetName} Score')
    main_ax.set_xlim(xmin,xmax)
    main_ax.set_ylim(ymin,ymax)
    main_ax.set_title(f'{Experiment}: {gsetName} Score in ATAC and GEX',pad=10)

    # Density plot for x below the main plot
    x_density_ax = fig.add_subplot(gs[2, 0], sharex=main_ax)
    sns.kdeplot(data=rna.obs[gsetName], color='#4daf4a', alpha=0.5, linewidth=2, ax=x_density_ax)
    x_density_ax.set_ylabel('Cell density')
    x_density_ax.set_xlabel(f'')
    x_density_ax.yaxis.tick_right()
    x_density_ax.yaxis.set_label_position("right")

    # Density plot for y left of the main plot
    y_density_ax = fig.add_subplot(gs[1, 1], sharey=main_ax)
	sns.kdeplot(data=atac.obs, color='#e41a1c', alpha=0.5, linewidth=2, ax=y_density_ax, y=gsetName)
    y_density_ax.set_ylabel('')
    y_density_ax.set_xlabel('Cell density')
    y_density_ax.xaxis.tick_bottom()

    # Hide x and y ticks for density plots
    plt.setp(x_density_ax.get_xticklabels(), visible=False)
    plt.setp(y_density_ax.get_yticklabels(), visible=False)

    # Add color bar to the main plot
    #cax = fig.add_subplot(gs[1, 2])
    #cbar = fig.colorbar(main_ax.contourf(x, y, z, cmap='plasma', levels=30, antialiased=True), cax=cax)
    #cbar.set_label('Density', rotation=90)

    #plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.tight_layout()
    if dens:
    	plt.savefig(f'{Experiment}_GEX_ATAC_{gsetName}_Score_Density.pdf')
    else:
    	plt.savefig(f'{Experiment}_GEX_ATAC_{gsetName}_Score_Scatter.pdf')
    plt.close()

for gset_name in gene_sigs:
	print(f'{Experiment}: {gset_name} scatter plot')
	scatter_gene_sigs(rna,atac,gset_name,dens=False)
	print(f'{Experiment}: {gset_name} density plot')
	scatter_gene_sigs(rna,atac,gset_name,dens=True)

###
### Cluster analysis pan samples
###
type_col_name='TypeDiffAmbig'

print(f'Reading {Experiment} H5ad files...')
#rna = sc.read_h5ad(DirRNA+ Exp_RNA + annot_ext)
#atac = sc.read_h5ad(DirATAC+ Exp_ATAC + annot_ext)

GEXdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': rna.obs_names, 'type': rna.obs[type_col_name]})
GEXdf['type'] = GEXdf['type'].str.split().str[0]
ATACdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': atac.obs_names, 'type': atac.obs[type_col_name]})
ATACdf['type'] = ATACdf['type'].str.split().str[0]

merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')

#commonBClist = list(set(atac.obs.index).intersection(set(rna.obs.index)))
#print(np.all(np.array([True for el in commonBClist if el in list(merged_df['obs_names'])])),'Barcode Intersection')
#print(f'RNA lost {len(GEXdf)-len(merged_df)} cells\nATAC lost {len(ATACdf)-len(merged_df)}')
print(f'Plotting {Experiment}')

#########
#Matrix comparison
#########

def jaccard_types(merged_df,subset_name,Experiment=Experiment):
	groups1 = np.array(merged_df['type_x'])
	cl1 = np.unique(groups1)
	groups2 = np.array(merged_df['type_y'])
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
	ax.set_xticklabels(cl2[sorted_indices_x],rotation=45,rotation_mode='anchor',ha='right')
	ax.set_yticklabels(cl1[sorted_indices_y])

	# Add the number of observations in each group on x and y axes
	for i, (value, label) in enumerate(zip(cl2[sorted_indices_x], ax.get_xticklabels())):
	    ax.text(i, -0.56, f'{len(np.where(groups2 == value)[0])}', ha='left', va='center',rotation=45,rotation_mode='anchor')

	for i, (value, label) in enumerate(zip(cl1[sorted_indices_y], ax.get_yticklabels())):
	    ax.text(-0.45, i, f'{len(np.where(groups1 == value)[0])}', ha='left', va='center')
	    label.set_y(i - 0.3)  # Adjust the position of the current label

	ax.set_xlabel('GEX Cell Types')
	ax.set_ylabel('ATAC Cell Types')
	if len(subset_name)>0:
		ax.set_title(f'{subset_name[:-3]}',pad=30)
	else:
		ax.set_title(f'{Experiment}',pad=30)
	divider = make_axes_locatable(plt.gca())
	cax = divider.append_axes("right", "5%", pad="3%")
	plt.colorbar(im, cax=cax)
	plt.tight_layout()

	# save the plot
	plt.savefig(f'./figures/{Experiment}_{subset_name}_Jaccard.pdf')
	plt.close()


#########
#Matrix comparison Pan Sample
#########

jaccard_types(merged_df,'ALL',Experiment=Experiment)


#########
#Matrix comparison Per Sample
#########
for sample in rna.obs['batch'].cat.categories:
	print(f'Jaccard per sample ({sample})')
	dataGEX = rna[rna.obs['batch']==sample].copy()
	dataATAC = atac[atac.obs['batch']==sample].copy()

	GEXdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': dataGEX.obs_names, 'type': dataGEX.obs[type_col_name]})
	GEXdf['type'] = GEXdf['type'].str.split().str[0]
	ATACdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': dataATAC.obs_names, 'type': dataATAC.obs[type_col_name]})
	ATACdf['type'] = ATACdf['type'].str.split().str[0]
	merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')
	jaccard_types(merged_df,f'{sample}_SP')


#########
#Matrix comparison Per timepoints
#########

for timepoint in rna.obs['timePoint'].cat.categories:
	print(f'Jaccard per timepoint ({timepoint})')
	dataGEX = rna[rna.obs['timePoint']==timepoint].copy()
	dataATAC = atac[atac.obs['timePoint']==timepoint].copy()

	GEXdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': dataGEX.obs_names, 'type': dataGEX.obs[type_col_name]})
	GEXdf['type'] = GEXdf['type'].str.split().str[0]
	ATACdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': dataATAC.obs_names, 'type': dataATAC.obs[type_col_name]})
	ATACdf['type'] = ATACdf['type'].str.split().str[0]
	merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')
	jaccard_types(merged_df,f'{timepoint}_TP')


# Concat samples
pdfjam --nup 4x3 --landscape --a4paper *ENZ*SP_Jaccard.pdf --outfile Wouter21_ENZ_SamplesJaccard.pdf
pdfjam --nup 1x2 --landscape --a4paper *ENZ*TP_Jaccard.pdf --outfile Wouter21_ENZ_TimepointJaccard.pdf
rm *ENZ*P_Jaccard.pdf


pdfjam --nup 4x3 --landscape --a4paper *SING*SP_Jaccard.pdf --outfile Wouter21_SING_SamplesJaccard.pdf
pdfjam --nup 3x2 --landscape --a4paper *SING*TP_Jaccard.pdf --outfile Wouter21_SING_TimepointJaccard.pdf
rm *SING*P_Jaccard.pdf



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



################################################################################################
################################################################################################
################################################################################################
# Check new sigs 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde
import os
import plotly.graph_objects as go

### TO CHECK WHEN CHANGING SAMPLES ###
met = 'SING'#,'SING']:
### TO CHECK WHEN CHANGING SAMPLES ###


refDir = '/mnt/etemp/ahrmad/wouter/refs'
DirRNA = f'/mnt/etemp/ahrmad/wouter/RNA_CB/annot_other_norm/'
DirATAC = f'/mnt/etemp/ahrmad/wouter/RNA_CB/annot_other_norm/'
Exp_ATAC = f'Wouter21_{met}_CB_annot_All_newSigs.h5ad'
Exp_RNA = f'Wouter21_{met}_CB_annot_All.h5ad'
Experiment=f'Wouter21_{met}'
annot_ext = ''
gene_set_table = f'{refDir}/table_s8_summary.txt'
gene_sigs = ['Epi_Luminal_1', 'Epi_Luminal_2Psca']

assert (os.path.exists(DirRNA+ Exp_RNA + annot_ext)) & (os.path.exists(DirATAC+ Exp_ATAC + annot_ext)),'ScoreGene H5ADs missing'

rna = sc.read_h5ad(DirRNA+ Exp_RNA + annot_ext)
atac = sc.read_h5ad(DirATAC+ Exp_ATAC + annot_ext)


# Make sure the same cells are in both adatas
print(f'keeping cells in both modalities')
rna.obs = rna.obs.reindex(index=sorted(rna.obs.index))
atac.obs = atac.obs.reindex(index=sorted(atac.obs.index))
lenATAC,lenGEX = atac.n_obs,rna.n_obs
atac = atac[atac.obs.index.isin(rna.obs.index)]
rna = rna[rna.obs.index.isin(atac.obs.index)]
assert (rna.obs.index == atac.obs.index).all()
print(f'rna lost {lenGEX-rna.n_obs} cells\natac lost {lenATAC-atac.n_obs} cells')


type_col_name='Annotation'

print(f'Reading {Experiment} H5ad files...')
#rna = sc.read_h5ad(DirRNA+ Exp_RNA + annot_ext)
#atac = sc.read_h5ad(DirATAC+ Exp_ATAC + annot_ext)

GEXdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': rna.obs_names, 'type': rna.obs[type_col_name]})
GEXdf['type'] = GEXdf['type'].str.split().str[0]
ATACdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': atac.obs_names, 'type': atac.obs[type_col_name]})
ATACdf['type'] = ATACdf['type'].str.split().str[0]

merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')

jaccard_types(merged_df,'ALL',Experiment=Experiment)

for timepoint in rna.obs['timePoint'].cat.categories:
	print(f'Jaccard per timepoint ({timepoint})')
	dataGEX = rna[rna.obs['timePoint']==timepoint].copy()
	dataATAC = atac[atac.obs['timePoint']==timepoint].copy()

	GEXdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': dataGEX.obs_names, 'type': dataGEX.obs[type_col_name]})
	GEXdf['type'] = GEXdf['type'].str.split().str[0]
	ATACdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': dataATAC.obs_names, 'type': dataATAC.obs[type_col_name]})
	ATACdf['type'] = ATACdf['type'].str.split().str[0]
	merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')
	jaccard_types(merged_df,f'{timepoint}_TP')



